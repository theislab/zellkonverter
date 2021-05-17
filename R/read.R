#' Read H5AD
#'
#' Reads a H5AD file and returns a \linkS4class{SingleCellExperiment} object.
#'
#' @param file String containing a path to a `.h5ad` file.
#' @param X_name Name used when saving `X` as an assay. If `NULL` looks for an
#' `X_name` value in `uns`, otherwise uses `"X"`.
#' @param use_hdf5 Logical scalar indicating whether assays should be
#' loaded as HDF5-based matrices from the **HDF5Array** package.
#'
#' @details
#' Setting `use_hdf5 = TRUE` allows for very large datasets to be efficiently
#' represented on machines with little memory. However, this comes at the cost
#' of access speed as data needs to be fetched from the HDF5 file upon request.
#'
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies. This will not be performed on
#' any subsequent run or if any other **zellkonverter** function has been run
#' prior to this one.
#'
#' @return A \linkS4class{SingleCellExperiment} object is returned.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
#' sce <- readH5AD(file)
#' class(assay(sce))
#'
#' sce2 <- readH5AD(file, use_hdf5 = TRUE)
#' class(assay(sce2))
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @seealso
#' [`writeH5AD()`], to write a \linkS4class{SingleCellExperiment} object to a
#' H5AD file.
#'
#' [`AnnData2SCE()`], for developers to convert existing AnnData instances to a
#' \linkS4class{SingleCellExperiment}.
#'
#' @export
#' @importFrom basilisk basiliskRun
readH5AD <- function(file, X_name = NULL, use_hdf5 = FALSE) {
    file <- path.expand(file)

    basiliskRun(
        env = zellkonverterAnnDataEnv,
        fun = .H5ADreader,
        file = file,
        X_name = X_name,
        backed = use_hdf5
    )
}

#' @importFrom reticulate import
.H5ADreader <- function(file, X_name = NULL, backed = FALSE) {
    anndata <- import("anndata")
    adata <- anndata$read_h5ad(file, backed = if (backed) "r" else FALSE)
    AnnData2SCE(adata, X_name = X_name, hdf5_backed = backed)
}

#' @importFrom S4Vectors I DataFrame wmsg
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<- colPairs<- rowPairs<-
.native_reader <- function(file, backed=FALSE) {
    contents <- .list_contents(file)

    # Let's read in the X matrix first... if it's there.
    if (!"X" %in% names(contents)) {
        stop("missing an 'X' entry in '", file, "'")
    }
    all.assays <- list(X = .read_matrix(file, "X", contents[["X"]], backed=backed))

    for (l in names(contents[["layers"]])) {
        tryCatch({
            all.assays[[l]] <- .read_matrix(file, file.path("layers", l), contents[["layers"]][[l]], backed=backed)
        }, error=function(e) {
            warning(wmsg("setting additional assays from 'layers' failed for '", file, "':\n  ", conditionMessage(e)))
        })
    }

    sce <- SingleCellExperiment(all.assays)

    # Adding the various pieces of data.
    tryCatch({
        cd <- .read_dim_data(file, "obs", contents[["obs"]])
        if (!is.null(cd)) {
            colData(sce) <- cd
        }
    }, error=function(e) {
        warning(wmsg("setting 'colData' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    tryCatch({
        rd <- .read_dim_data(file, "var", contents[["var"]])
        if (!is.null(rd)) {
            rowData(sce) <- rd
        }
    }, error=function(e) {
        warning(wmsg("setting 'rowData' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    # Adding the reduced dimensions and other bits and pieces.
    tryCatch({
        reducedDims(sce) <- .read_dim_mats(file, "obsm", contents[["obsm"]]) 
    }, error=function(e) {
        warning(wmsg("setting 'reducedDims' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    tryCatch({
        row.mat <- .read_dim_mats(file, "varm", contents[["varm"]])
        if (length(row.mat)) {
            row.mat.df <- do.call(DataFrame, lapply(row.mat, I))
            rowData(sce) <- cbind(rowData(sce), row.mat.df)
        }
    }, error=function(e) {
        warning(wmsg("extracting 'varm' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    # Adding pairings, if any exist.
    tryCatch({
        rowPairs(sce) <- .read_dim_pairs(file, "varp", contents[["varp"]])
    }, error=function(e) {
        warning(wmsg("setting 'rowPairs' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    tryCatch({
        colPairs(sce) <- .read_dim_pairs(file, "obsp", contents[["obsp"]])
    }, error=function(e) {
        warning(wmsg("setting 'colPairs' failed for '", file, "':\n  ", conditionMessage(e)))
    })

    if ("uns" %in% names(contents)) {
        tryCatch({
            metadata(sce) <- rhdf5::h5read(file, "uns")
        }, error=function(e) {
            warning(wmsg("setting 'metadata' failed for '", file, "':\n  ", conditionMessage(e)))
        })
    }

    sce
}

.list_contents <- function(file) {
    manifest <- rhdf5::h5ls(file)

    set_myself <- function(x, series, value) {
        if (length(series)!=1) {
            value <- set_myself(x[[series[1]]], series[-1], value)
        }
        if (is.null(x)) {
            x <- list()
        }
        x[[series[1]]] <- value
        x
    }

    contents <- list()
    for (i in seq_len(nrow(manifest))) {
        components <- c(strsplit(manifest[i, "group"], "/")[[1]], manifest[i, "name"])
        if (components[1] == "") {
            components <- components[-1]
        }

        info <- manifest[i, c('otype', 'dclass', 'dim')]
        if (info$otype=="H5I_GROUP") {
            info <- list()
        }
        contents <- set_myself(contents, components, info)
    }

    contents
}

.read_matrix <- function(file, path, fields, backed) {
    if (is.data.frame(fields)) {
        mat <- HDF5Array::HDF5Array(file, path)
    } else {
        mat <- HDF5Array::H5SparseMatrix(file, path)
    }
    if (!backed) {
         if (DelayedArray::is_sparse(mat)) {
             mat <- as(mat, "sparseMatrix")
         } else {
             mat <- as.matrix(mat)
         }
    }
    mat
}

#' @importFrom S4Vectors DataFrame
.read_dim_data <- function(file, path, fields) {
    col.names <- setdiff(names(fields), c("__categories", "_index"))
    out.cols <- list()
    for (i in col.names) {
        out.cols[[i]] <- as.vector(rhdf5::h5read(file, file.path(path, i)))
    }

    cat.names <- names(fields[["__categories"]])
    for (i in cat.names) {
        levels <- as.vector(rhdf5::h5read(file, file.path(path, "__categories", i)))
        out.cols[[i]] <- factor(out.cols[[i]], levels)
    }

    if (!is.null(fields[["_index"]])) {
        indices <- as.vector(rhdf5::h5read(file, file.path(path, "_index")))
    } else {
        indices <- NULL
    }

    if (length(out.cols)) {
        df <- do.call(DataFrame, out.cols)
        rownames(df) <- indices
    } else if (!is.null(indices)) {
        df <- DataFrame(row.names=indices)
    } else {
        df <- NULL
    }

    df
}

.read_dim_mats <- function(file, path, fields) {
    all.contents <- list()
    for (i in names(fields)) {
        # because everything's transposed.
        all.contents[[i]] <- t(rhdf5::h5read(file, file.path(path, i)))
    }
    all.contents
}

.read_dim_pairs <- function(file, path, fields) {
    all.pairs <- list()
    for (i in names(fields)) {
        mat <- HDF5Array::H5SparseMatrix(file, file.path(path, i))
        all.pairs[[i]] <- as(mat, "sparseMatrix")
    }
    all.pairs
}
