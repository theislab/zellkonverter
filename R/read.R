#' Read H5AD
#'
#' Reads a H5AD file and returns a \linkS4class{SingleCellExperiment} object.
#'
#' @param file String containing a path to a `.h5ad` file.
#' @param X_name Name used when saving `X` as an assay. If `NULL` looks for an
#' `X_name` value in `uns`, otherwise uses `"X"`.
#' @param use_hdf5 Logical scalar indicating whether assays should be
#' loaded as HDF5-based matrices from the **HDF5Array** package.
#' @param reader Which HDF5 reader to use. Either `"python"` for reading with
#' the **anndata** Python package via **reticulate** or `"R"` for
#' **zellkonverter**'s native R reader.
#' @param version A string giving the version of the **anndata** Python library
#' to use. Allowed values are available in `.AnnDataVersions`. By default the
#' latest version is used.
#' @param verbose Logical scalar indicating whether to print progress messages.
#' If `NULL` uses `getOption("zellkonverter.verbose")`.
#' @inheritDotParams AnnData2SCE -adata -hdf5_backed
#'
#' @details
#' Setting `use_hdf5 = TRUE` allows for very large datasets to be efficiently
#' represented on machines with little memory. However, this comes at the cost
#' of access speed as data needs to be fetched from the HDF5 file upon request.
#'
#' Setting `reader = "R"` will use an experimental native R reader instead of
#' reading the file into Python and converting the result. This avoids the need
#' for a Python environment and some of the issues with conversion but is still
#' under development and is likely to return slightly different output.
#'
#' See [AnnData-Environment] for more details on **zellkonverter** Python
#' environments.
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
#'
#' sce3 <- readH5AD(file, reader = "R")
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
#' @importFrom methods slot
readH5AD <- function(file, X_name = NULL, use_hdf5 = FALSE,
                     reader = c("python", "R"), version = NULL, verbose = NULL,
                     ...) {
    file <- path.expand(file)
    reader <- match.arg(reader)

    if (reader == "python") {
        .ui_info("Using the {.field Python} reader")
        env <- zellkonverterAnnDataEnv(version)
        version <- gsub("zellkonverterAnnDataEnv-", "", slot(env, "envname"))
        .ui_info("Using {.field anndata} version {.field {version}}")

        sce <- basiliskRun(
            env = env,
            fun = .H5ADreader,
            file = file,
            X_name = X_name,
            backed = use_hdf5,
            verbose = verbose,
            ...
        )

    } else if (reader == "R") {
        sce <- .native_reader(file, backed = use_hdf5, verbose = verbose)
    }

    return(sce)
}

#' @importFrom reticulate import
.H5ADreader <- function(file, X_name = NULL, backed = FALSE, verbose = NULL, ...) {
    anndata <- import("anndata")
    .ui_step(
        "Reading {.file { .trim_path(file)} }",
        msg_done = "Read {.file { .trim_path(file) }}",
        spinner = TRUE
    )
    adata <- anndata$read_h5ad(file, backed = if (backed) "r" else FALSE)
    cli::cli_progress_done()
    AnnData2SCE(adata, X_name = X_name, hdf5_backed = backed, verbose = verbose,
                ...)
}

#' @importFrom S4Vectors I DataFrame wmsg
#' @importFrom SummarizedExperiment assays assays<- rowData colData rowData<- colData<-
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<- colPairs<- rowPairs<-
.native_reader <- function(file, backed = FALSE, verbose = FALSE) {
    .ui_info("Using the {.field R} reader")
    .ui_step("Reading {.file {file}}", spinner = TRUE)

    contents <- .list_contents(file)

    all.assays <- list()

    # Let's read in the X matrix first... if it's there.
    if ("X" %in% names(contents)) {
        all.assays[["X"]] <- .read_matrix(file, "X", contents[["X"]], backed = backed)
    }

    for (layer in names(contents[["layers"]])) {
        tryCatch(
            {
                all.assays[[layer]] <- .read_matrix(
                    file,
                    file.path("layers", layer),
                    contents[["layers"]][[layer]],
                    backed = backed
                )
            },
            error = function(e) {
                warning(wmsg(
                    "setting additional assays from 'layers' failed for '",
                    file, "':\n  ", conditionMessage(e)
                ))
            }
        )
    }

    sce <- SingleCellExperiment(all.assays)

    # Adding the various pieces of data.
    tryCatch(
        {
            col_data <- .read_dim_data(file, "obs", contents[["obs"]])
            if (!is.null(col_data)) {
                colData(sce) <- col_data
            }
        },
        error = function(e) {
            warning(wmsg(
                "setting 'colData' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    tryCatch(
        {
            row_data <- .read_dim_data(file, "var", contents[["var"]])
            if (!is.null(row_data)) {
                rowData(sce) <- row_data
                # Manually set SCE rownames, because setting rowData
                # doesn't seem to set them. (Even tho setting colData
                # does set the colnames)
                rownames(sce) <- rownames(row_data)
            }
        },
        error = function(e) {
            warning(wmsg(
                "setting 'rowData' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    # Adding the reduced dimensions and other bits and pieces.
    tryCatch(
        {
            reducedDims(sce) <- .read_dim_mats(file, "obsm", contents[["obsm"]])
        },
        error = function(e) {
            warning(wmsg(
                "setting 'reducedDims' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    tryCatch(
        {
            row_mat <- .read_dim_mats(file, "varm", contents[["varm"]])
            if (length(row_mat)) {
                row_mat_df <- do.call(DataFrame, lapply(row_mat, I))
                rowData(sce) <- cbind(rowData(sce), row_mat_df)
            }
        },
        error = function(e) {
            warning(wmsg(
                "extracting 'varm' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    # Adding pairings, if any exist.
    tryCatch(
        {
            rowPairs(sce) <- .read_dim_pairs(file, "varp", contents[["varp"]])
        },
        error = function(e) {
            warning(wmsg(
                "setting 'rowPairs' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    tryCatch(
        {
            colPairs(sce) <- .read_dim_pairs(file, "obsp", contents[["obsp"]])
        },
        error = function(e) {
            warning(wmsg(
                "setting 'colPairs' failed for '", file, "':\n  ",
                conditionMessage(e)
            ))
        }
    )

    if ("uns" %in% names(contents)) {
        tryCatch(
            {
                uns <- rhdf5::h5read(file, "uns")
                uns <- .convert_element(
                    uns, "uns", file, recursive=TRUE
                )
                metadata(sce) <- uns
            },
            error = function(e) {
                warning(wmsg(
                    "setting 'metadata' failed for '", file, "':\n  ",
                    conditionMessage(e)
                ))
            }
        )
    }

    if (("X_name" %in% names(metadata(sce))) && ("X" %in% names(contents))) {
        stopifnot(names(assays(sce))[1] == "X") #should be true b/c X is read 1st
        names(assays(sce))[1] <- metadata(sce)[["X_name"]]
        metadata(sce)[["X_name"]] <- NULL
    }

    sce
}

.list_contents <- function(file) {
    manifest <- rhdf5::h5ls(file)

    set_myself <- function(x, series, value) {
        if (length(series) != 1) {
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
        components <- c(
            strsplit(manifest[i, "group"], "/")[[1]],
            manifest[i, "name"]
        )
        if (components[1] == "") {
            components <- components[-1]
        }

        info <- manifest[i, c("otype", "dclass", "dim")]
        if (info$otype == "H5I_GROUP") {
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

.convert_element <- function(obj, path, file, recursive=FALSE) {
    element_attrs <- rhdf5::h5readAttributes(file, path)

    # Convert categorical element for AnnData v0.8+
    if (identical(element_attrs[["encoding-type"]], "categorical") &&
        all(c("codes", "categories") %in% names(obj))) {
        codes <- obj[["codes"]] + 1
        codes[codes == 0] <- NA
        levels <- obj[["categories"]]

        ord <- as.logical(element_attrs[["ordered"]])

        obj <- factor(levels[codes], levels=levels, ordered=ord)
        return(obj)
    }

    # Handle booleans. Non-nullable booleans have encoding-type
    # "array", so we have to infer the type from the enum levels
    if (is.factor(obj) && identical(levels(obj), c("FALSE", "TRUE"))) {
        obj <- as.logical(obj)
        return(obj)
    }

    # Recursively convert element members
    if (recursive && is.list(obj) && !is.null(names(obj))) {
        for (k in names(obj)) {
            obj[[k]] <- rhdf5::h5read(file, file.path(path, k))
            obj[[k]] <- .convert_element(
                obj[[k]], file.path(path, k),
                file, recursive=TRUE
            )
        }
    }

    if (is.list(obj) && !is.null(names(obj))) {
        names(obj) <- make.names(names(obj))
    }

    obj
}

#' @importFrom S4Vectors DataFrame
.read_dim_data <- function(file, path, fields) {
    col_names <- setdiff(names(fields), "__categories")
    out_cols <- list()
    for (col_name in col_names) {
        vec <- rhdf5::h5read(file, file.path(path, col_name))

        vec <- .convert_element(
            vec, file.path(path, col_name),
            file, recursive=FALSE
        )

        if (!is.factor(vec)) {
            vec <- as.vector(vec)
        }

        out_cols[[col_name]] <- vec
    }

    # for AnnData versions <= 0.7
    cat_names <- names(fields[["__categories"]])
    for (cat_name in cat_names) {
        levels <- as.vector(
            rhdf5::h5read(file, file.path(path, "__categories", cat_name))
        )
        out_cols[[cat_name]] <- factor(out_cols[[cat_name]])
        levels(out_cols[[cat_name]]) <- levels
    }

    ## rhdf5::h5readAttributes(file, "var") |> str()
    ## List of 4
    ##  $ _index          : chr "feature_id"
    ##  $ column-order    : chr [1:4(1d)] "feature_is_filtered" "feature_name" "feature_reference" "feature_biotype"
    ##  $ encoding-type   : chr "dataframe"
    ##  $ encoding-version: chr "0.2.0"
    attributes <- rhdf5::h5readAttributes(file, path)
    index <- attributes[["_index"]]
    if (!is.null(index)) {
        indices <- out_cols[[index]]
    } else {
        indices <- NULL
    }

    column_order <- attributes[["column-order"]]
    if (!is.null(column_order)) {
        out_cols <- out_cols[column_order]
    }

    if (length(out_cols)) {
        df <- do.call(DataFrame, out_cols)
        rownames(df) <- indices
    } else if (!is.null(indices)) {
        df <- DataFrame(row.names = indices)
    } else {
        df <- NULL
    }

    df
}

.read_dim_mats <- function(file, path, fields) {
    all.contents <- list()
    for (field in names(fields)) {
        # because everything's transposed.
        all.contents[[field]] <- t(rhdf5::h5read(file, file.path(path, field)))
    }
    all.contents
}

.read_dim_pairs <- function(file, path, fields) {
    all.pairs <- list()
    for (field in names(fields)) {
        mat <- HDF5Array::H5SparseMatrix(file, file.path(path, field))
        all.pairs[[field]] <- as(mat, "sparseMatrix")
    }
    all.pairs
}
