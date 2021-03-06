#' Convert AnnData between and SingleCellExperiment
#'
#' Conversion between Python AnnData objects and
#' \linkS4class{SingleCellExperiment} objects.
#'
#' @details
#' These functions assume that an appropriate Python environment has already
#' been loaded. As such, they are largely intended for developer use, most
#' typically inside a **basilisk** context.
#'
#' The conversion is not entirely lossless. The current mapping is shown below
#' (also at <https://tinyurl.com/AnnData2SCE>):
#'
#' \if{html}{
#'     \figure{AnnData2SCE.png}{options: width=800, alt="SCE-AnnData map"}
#' }
#' \if{latex}{\figure{AnnData2SCE.png}{options: width=5in}}
#'
#' In `SCE2AnnData()`, matrices are converted to a **numpy**-friendly format.
#' Sparse matrices are converted to \linkS4class{dgCMatrix} objects while all
#' other matrices are converted into ordinary matrices. If `skip_assays = TRUE`,
#' empty sparse matrices are created instead and the user is expected to fill in
#' the assays on the Python side.
#'
#' For `AnnData2SCE()`, a warning is raised if there is no corresponding R
#' format for a matrix in the AnnData object, and an empty sparse matrix is
#' created instead as a placeholder. If `skip_assays = NA`, no warning is
#' emitted but variables are created in the [`int_metadata()`] of the output to
#' specify which assays were skipped.
#'
#' If `skip_assays = TRUE`, empty sparse matrices are created for all assays,
#' regardless of whether they might be convertible to an R format or not.
#' In both cases, the user is expected to fill in the assays on the R side,
#' see [`readH5AD()`] for an example.
#'
#' We attempt to convert between items in the \linkS4class{SingleCellExperiment}
#' [`metadata()`] slot and the `AnnData` `uns` slot. If an item cannot be
#' converted a warning will be raised.
#'
#' Values stored in the `varm` slot of an `AnnData` object are stored in a
#' column of [`rowData()`] in a \linkS4class{SingleCellExperiment}
#' as a \linkS4class{DataFrame} of matrices. If this column is present an
#' attempt is made to transfer this information when converting from
#' \linkS4class{SingleCellExperiment} to `AnnData`.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @return `AnnData2SCE()` will return a \linkS4class{SingleCellExperiment}
#' containing the equivalent data from `adata`.
#'
#' `SCE2AnnData()` will return a **reticulate** reference to an AnnData object
#' containing the content of `sce`.
#'
#' @seealso
#' [`writeH5AD()`] and [`readH5AD()`] for dealing directly with H5AD files.
#'
#' @examples
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(basilisk)
#'     library(scRNAseq)
#'     seger <- SegerstolpePancreasData()
#'
#'     # These functions are designed to be run inside
#'     # a specified Python environment
#'     roundtrip <- basiliskRun(fun = function(sce) {
#'         # Convert SCE to AnnData:
#'         adata <- zellkonverter::SCE2AnnData(sce)
#'
#'         # Maybe do some work in Python on 'adata':
#'         # BLAH BLAH BLAH
#'
#'         # Convert back to an SCE:
#'         zellkonverter::AnnData2SCE(adata)
#'     }, env = zellkonverterAnnDataEnv, sce = seger)
#' }
#' @name AnnData-Conversion
#' @rdname AnnData-Conversion
NULL

#' @rdname AnnData-Conversion
#'
#' @param adata A **reticulate** reference to a Python AnnData object.
#' @param skip_assays Logical scalar indicating whether to skip conversion of
#' any assays in `sce` or `adata`, replacing them with empty sparse matrices
#' instead.
#' @param hdf5_backed Logical scalar indicating whether HDF5-backed matrices
#' in `adata` should be represented as HDF5Array objects. This assumes that
#' `adata` is created with `backed="r"`.
#'
#' @export
#' @importFrom methods selectMethod is
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom reticulate import_builtins
AnnData2SCE <- function(adata, X_name = NULL, skip_assays = FALSE,
                        hdf5_backed = TRUE) {
    py_builtins <- import_builtins()

    dims <- unlist(adata$shape)
    dims <- rev(dims)

    x_out <- .extract_or_skip_assay(
        skip_assays = skip_assays,
        hdf5_backed = hdf5_backed,
        dims = dims,
        mat = adata$X,
        name = "'X' matrix"
    )

    meta_list <- .convert_anndata_slot(
        adata, "uns", py_builtins$list(adata$uns$keys())
    )

    x_mat <- x_out$mat
    colnames(x_mat) <- adata$obs_names$to_list()
    rownames(x_mat) <- adata$var_names$to_list()
    skipped_x <- x_out$skipped

    if (is.null(X_name)) {
        if ("X_name" %in% names(meta_list)) {
            X_name <- meta_list[["X_name"]]
            message("Note: Using stored X_name value '", X_name, "'")
            meta_list[["X_name"]] <- NULL
        } else {
            X_name <- "X"
        }
    }

    assays_list <- list()
    assays_list[[X_name]] <- x_mat

    layer_names <- names(py_builtins$dict(adata$layers))
    skipped_layers <- character(0)
    for (layer_name in layer_names) {
        layer_out <- .extract_or_skip_assay(
            skip_assays = skip_assays,
            hdf5_backed = hdf5_backed,
            dims = dims,
            mat = adata$layers$get(layer_name),
            name = sprintf("'%s' layer matrix", layer_name)
        )
        if (layer_out$skipped) {
            skipped_layers <- c(skipped_layers, layer_name)
        }
        assays_list[[layer_name]] <- layer_out$mat
    }

    varp_list <- .convert_anndata_slot(
        adata, "varp", py_builtins$list(adata$varp$keys())
    )

    obsp_list <- .convert_anndata_slot(
        adata, "obsp", py_builtins$list(adata$obsp$keys())
    )

    row_data <- DataFrame(adata$var)

    varm_list <- .convert_anndata_slot(adata, "varm", adata$varm_keys())

    if (length(varm_list) > 0) {
        # Create an empty DataFrame with the correct number of rows
        varm_df <- make_zero_col_DFrame(adata$n_vars)
        for (varm_name in names(varm_list)) {
            varm_df[[varm_name]] <- varm_list[[varm_name]]
        }
        row_data$varm <- varm_df
    }

    output <- SingleCellExperiment(
        assays = assays_list,
        rowData = row_data,
        colData = adata$obs,
        reducedDims = py_builtins$dict(adata$obsm),
        metadata = meta_list,
        rowPairs = varp_list,
        colPairs = obsp_list
    )

    # Specifying which assays got skipped, if the skipping was variable.
    if (is.na(skip_assays)) {
        int_metadata(output)$skipped_x <- skipped_x
        int_metadata(output)$skipped_layers <- skipped_layers
    }

    if (length(varm_list) > 0) {
        int_metadata(output)$has_varm <- names(varm_list)
    }

    output
}

#' @importFrom Matrix t
.extract_or_skip_assay <- function(skip_assays, hdf5_backed, dims, mat, name) {
    skipped <- FALSE

    if (isTRUE(skip_assays)) {
        # Value of 'mat' is never used so the promise never evaluates; thus,
        # skip_assays=TRUE avoids any actual transfer of content from Python.
        mat <- .make_fake_mat(dims)
    } else {
        if (hdf5_backed && is(mat, "python.builtin.object")) {
            file <- as.character(mat$file$id$name)
            name <- as.character(mat$name)
            if (.h5isgroup(file, name)) {
                mat <- HDF5Array::H5SparseMatrix(file, name)
            } else {
                mat <- HDF5Array::HDF5Array(file, name)
            }
        } else {
            mat <- try(t(mat), silent = TRUE)
            if (is(mat, "try-error")) {
                if (isFALSE(skip_assays)) {
                    warning(
                        name,
                        " does not support transposition and has been skipped"
                    )
                }
                mat <- .make_fake_mat(dims)
                skipped <- TRUE
            }
        }
    }

    if (is(mat, "dgRMatrix")) {
        mat <- as(mat, "CsparseMatrix")
    }

    list(mat = mat, skipped = skipped)
}

#' @importFrom Matrix sparseMatrix
.make_fake_mat <- function(dims) {
    sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = dims
    )
}

# Borrowed from HDF4Array
# https://github.com/Bioconductor/HDF5Array/blob/fb015cf1c789bbb905a6bf8af2c7b50a24a60795/R/h5utils.R#L66-L75
.h5isgroup <- function(filepath, name) {
    fid <- rhdf5::H5Fopen(filepath, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(fid))
    gid <- try(rhdf5::H5Gopen(fid, name), silent = TRUE)
    ans <- !inherits(gid, "try-error")
    if (ans) {
        rhdf5::H5Gclose(gid)
    }

    return(ans)
}

.convert_anndata_slot <- function(adata, slot_name, slot_keys) {
    py_builtins <- import_builtins()

    converted_list <- list()

    for (key in slot_keys) {
        tryCatch(
            {
                item <- adata[slot_name][[key]]

                item_type <- py_builtins$str(py_builtins$type(item))
                if (grepl("OverloadedDict", item_type)) {
                    item <- py_builtins$dict(item)
                }

                if (is(item, "python.builtin.object")) {
                    item <- reticulate::py_to_r(item)
                }

                converted_list[[key]] <- item
            },
            error = function(err) {
                warning(
                    "conversion failed for the item '",
                    key, "' in '", slot_name, "' with ",
                    "the following error and has been skipped\n",
                    "Conversion error message: ", err,
                    call. = FALSE
                )
            }
        )
    }

    return(converted_list)
}

#' @rdname AnnData-Conversion
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name For `SCE2AnnData()` name of the assay to use as the primary
#' matrix (`X`) of the AnnData object. If `NULL`, the first assay of `sce` will
#' be used by default. For `AnnData2SCE()` name used when saving `X` as an
#' assay. If `NULL` looks for an `X_name` value in `uns`, otherwise uses `"X"`.
#'
#' @export
#' @importFrom utils capture.output
#' @importFrom S4Vectors metadata
#' @importFrom reticulate import r_to_py
SCE2AnnData <- function(sce, X_name = NULL, skip_assays = FALSE) {
    anndata <- import("anndata")

    if (is.null(X_name)) {
        if (length(assays(sce)) == 0) {
            stop("'sce' does not contain any assays")
        }
        X_name <- assayNames(sce)[1]
        message("Note: using the '", X_name, "' assay as the X matrix")
    }

    if (!skip_assays) {
        X <- assay(sce, X_name)
        X <- .makeNumpyFriendly(X)
    } else {
        X <- fake_mat <- .make_fake_mat(rev(dim(sce)))
    }
    adata <- anndata$AnnData(X = X)

    col_data <- colData(sce)
    is_atomic <- vapply(col_data, is.atomic, NA)
    if (any(!is_atomic)) {
        non_atomic_cols <- colnames(col_data)[!is_atomic]
        warning(
            "The following colData columns are not atomic and will be stored ",
            "in metadata(sce)$.colData before conversion: ",
            paste(non_atomic_cols, collapse = ", ")
        )

        if (".colData" %in% names(metadata(sce))) {
            meta_list <- metadata(sce)$.colData
        } else {
            meta_list <- list()
        }

        for (col in non_atomic_cols) {
            store_name <- make.names(c(col, names(meta_list)), unique = TRUE)[1]
            meta_list[[store_name]] <- col_data[[col]]
        }

        col_data[non_atomic_cols] <- NULL
        metadata(sce)$.colData <- meta_list
    }

    if (ncol(col_data) > 0) {
        # Manually construct the data.frame to avoid mangling column names
        obs <- do.call(
            data.frame,
            c(
                as.list(col_data),
                check.names      = FALSE,
                stringsAsFactors = FALSE
            )
        )
        adata$obs <- obs
    }

    row_data <- rowData(sce)
    if (!is.null(int_metadata(sce)$has_varm)) {
        varm <- as.list(row_data[["varm"]])
        row_data[["varm"]] <- NULL
    }

    is_atomic <- vapply(row_data, is.atomic, NA)
    if (any(!is_atomic)) {
        non_atomic_cols <- colnames(row_data)[!is_atomic]
        warning(
            "The following rowData columns are not atomic and will be stored ",
            "in metadata(sce)$.rowData before conversion: ",
            paste(non_atomic_cols, collapse = ", ")
        )

        if (".rowData" %in% names(metadata(sce))) {
            meta_list <- metadata(sce)$.rowData
        } else {
            meta_list <- list()
        }

        for (col in non_atomic_cols) {
            store_name <- make.names(c(col, names(meta_list)), unique = TRUE)[1]
            meta_list[[store_name]] <- row_data[[col]]
        }

        row_data[non_atomic_cols] <- NULL
        metadata(sce)$.rowData <- meta_list
    }

    if (ncol(row_data) > 0) {
        # Manually construct the data.frame to avoid mangling column names
        var <- do.call(
            data.frame,
            c(
                as.list(row_data),
                check.names      = FALSE,
                stringsAsFactors = FALSE
            )
        )
        adata$var <- var
    }

    assay_names <- assayNames(sce)
    assay_names <- assay_names[!assay_names == X_name]
    if (length(assay_names) > 0) {
        if (!skip_assays) {
            assays_list <- assays(sce, withDimnames = FALSE)
            assays_list <- lapply(assays_list[assay_names], .makeNumpyFriendly)
        } else {
            assays_list <- rep(list(fake_mat), length(assay_names))
            names(assays_list) <- assay_names
        }
        adata$layers <- assays_list
    }

    red_dims <- as.list(reducedDims(sce))
    red_dims <- lapply(red_dims, .makeNumpyFriendly, transpose = FALSE)
    adata$obsm <- red_dims

    meta_list <- metadata(sce)
    meta_list <- .addListNames(meta_list)
    uns_list <- list()
    for (item_name in names(meta_list)) {
        item <- meta_list[[item_name]]
        tryCatch(
            {
                # Try to convert the item using reticulate, skip if it fails
                # Capture the object output printed by reticulate
                capture.output(r_to_py(item))
                uns_list[[item_name]] <- item
            },
            error = function(err) {
                warning(
                    "the '", item_name, "' item in 'metadata' cannot be ",
                    "converted to a Python type and has been skipped"
                )
            }
        )
    }
    uns_list[["X_name"]] <- X_name

    adata$uns <- reticulate::dict(uns_list)

    adata$varp <- as.list(rowPairs(sce, asSparse = TRUE))
    adata$obsp <- as.list(colPairs(sce, asSparse = TRUE))

    if (!is.null(int_metadata(sce)$has_varm)) {
        adata$varm <- varm
    }

    if (!is.null(colnames(sce))) {
        adata$obs_names <- colnames(sce)
    }

    if (!is.null(rownames(sce))) {
        adata$var_names <- rownames(sce)
    }

    adata
}

#' @importFrom methods as is
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom DelayedArray is_sparse
#' @importFrom Matrix t
.makeNumpyFriendly <- function(x, transpose = TRUE) {
    if (transpose) {
        x <- t(x)
    }

    # Code from Charlotte Soneson in kevinrue/velociraptor.
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}

.addListNames <- function(x) {
    if (length(x) == 0) {
        return(x)
    }

    if (is.null(names(x))) {
        names(x) <- paste0("item", seq_along(x))
        return(x)
    }

    list_names <- names(x)
    is_empty <- list_names == ""
    list_names[is_empty] <- paste0("item", seq_along(x)[is_empty])
    list_names <- make.names(list_names, unique = TRUE)

    names(x) <- list_names

    return(x)
}
