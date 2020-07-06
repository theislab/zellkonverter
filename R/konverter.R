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
#' \figure{AnnData2SCE.png}{options: width=800}
#'
#' In `SCE2AnnData()`, matrices are converted to a **numpy**-friendly format.
#' Sparse matrices are converted to \linkS4class{dgCMatrix} objects while all
#' other matrices are converted into ordinary matrices. If `skip_assays = TRUE`,
#' empty sparse matrices are created instead and the user is expected to fill in
#' the assays on the Python side.
#'
#' For `AnnData2SCE()`, a warning is raised if there is no corresponding R format
#' for a matrix in the AnnData object, and an empty sparse matrix is created 
#' instead as a placeholder. If `skip_assays = NA`, no warning is emitted
#' but variables are created in the [`int_metadata()`] of the output to specify
#' which assays were skipped.
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
#' as a \linkS4class{DataFrame} of matrices. No attempt is made to transfer this
#' information when converting from \linkS4class{SingleCellExperiment} to
#' `AnnData`.
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
#' library(basilisk)
#' library(scRNAseq)
#' seger <- SegerstolpePancreasData()
#'
#' # These functions are designed to be run inside
#' # a specified Python environment
#' roundtrip <- basiliskRun(fun = function(sce) {
#'      # Convert SCE to AnnData:
#'      adata <- SCE2AnnData(sce)
#'
#'      # Maybe do some work in Python on 'adata':
#'      # BLAH BLAH BLAH
#'
#'      # Convert back to an SCE:
#'      AnnData2SCE(adata)
#' }, env = zellkonverter:::anndata_env, sce = seger)
#'
#' @name AnnData-Conversion
#' @rdname AnnData-Conversion
NULL

#' @rdname AnnData-Conversion
#'
#' @param adata A **reticulate** reference to a Python AnnData object.
#' @param skip_assays Logical scalar indicating whether to skip conversion of
#' any assays in `sce` or `adata`, replacing them with empty sparse matrices
#' instead.
#'
#' @export
#' @importFrom methods selectMethod is
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom reticulate import_builtins
AnnData2SCE <- function(adata, skip_assays = FALSE) {
    py_builtins <- import_builtins()

    dims <- unlist(adata$shape)
    dims <- rev(dims)

    x_out <- .extract_or_skip_assay(
        skip_assays = skip_assays, 
        dims        = dims, 
        mat         = adata$X, 
        name        = "'X' matrix"
    )

    x_mat <- x_out$mat
    colnames(x_mat) <- adata$obs_names$to_list()
    rownames(x_mat) <- adata$var_names$to_list()
    skipped_x <- x_out$skipped

    assays_list <- list(X = x_mat)
    layer_names <- names(py_builtins$dict(adata$layers))
    skipped_layers <- character(0)

    for (layer_name in layer_names) {
        layer_out <- .extract_or_skip_assay(
            skip_assays = skip_assays,
            dims        = dims, 
            mat         = adata$layers$get(layer_name), 
            name        = sprintf("'%s' layer matrix", layer_name)
        )
        if (layer_out$skipped) {
            skipped_layers <- c(skipped_layers, layer_name)
        }
        assays_list[[layer_name]] <- layer_out$mat
    }

    uns_data <- adata$uns$data

    meta_list <- list()
    for (item_name in names(uns_data)) {
        item <- uns_data[[item_name]]
        if (!is(item, "python.builtin.object")) {
            meta_list[[item_name]] <- item
        } else {
            warning("the '", item_name, "' item in 'uns' cannot be converted ",
                    "to an R object and has been skipped")
        }
    }

    varp_list <- lapply(py_builtins$dict(adata$varp), function(v) v$todense())
    obsp_list <- lapply(py_builtins$dict(adata$obsp), function(v) v$todense())

    row_data <- DataFrame(adata$var)
    varm_list <- py_builtins$dict(adata$varm)
    if (length(varm_list) > 0) {
        # Create an empty DataFrame with the correct number of rows
        varm_df <- make_zero_col_DFrame(adata$n_vars)
        for (varm_name in names(varm_list)) {
            varm_df[[varm_name]] <- varm_list[[varm_name]]
        }
        row_data$varm <- varm_df
    }

    output <- SingleCellExperiment(
        assays      = assays_list,
        rowData     = row_data,
        colData     = adata$obs,
        reducedDims = py_builtins$dict(adata$obsm),
        metadata    = meta_list,
        rowPairs    = varp_list,
        colPairs    = obsp_list
    )

    # Specifying which assays got skipped, if the skipping was variable.
    if (is.na(skip_assays)) {
        int_metadata(output)$skipped_x <- skipped_x
        int_metadata(output)$skipped_layers <- skipped_layers
    }

    output
}

#' @importFrom Matrix t
.extract_or_skip_assay <- function(skip_assays, dims, mat, name) {
    skipped <- FALSE

    if (isTRUE(skip_assays)) {
        # Value of 'mat' is never used so the promise never evaluates; thus,
        # skip_assays=TRUE avoids any actual transfer of content from Python.
        mat <- .make_fake_mat(dims)
    } else {
        mat <- try(t(mat), silent=TRUE)
        if (is(mat, "try-error")) {
            if (isFALSE(skip_assays)) {
                warning(name, " does not support transposition and has been skipped")
            }
            mat <- .make_fake_mat(dims)
            skipped <- TRUE
        }
    }

    list(mat=mat, skipped=skipped)
}

#' @importFrom Matrix sparseMatrix
.make_fake_mat <- function(dims) {
    sparseMatrix(
        i    = integer(0),
        j    = integer(0),
        x    = numeric(0),
        dims = dims
    )
}

#' @rdname AnnData-Conversion
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name Name of the assay to use as the primary matrix (`X`) of the
#' AnnData object. If `NULL`, the first assay of `sce` will be used by default.
#'
#' @export
#' @importFrom Matrix t
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
        X <- t(assay(sce, X_name))
        X <- .makeNumpyFriendly(X)
    } else {
        X <- fake_mat <- .make_fake_mat(rev(dim(sce)))
    }
    adata <- anndata$AnnData(X = X)

    col_data <- colData(sce)
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
            assays_list <- lapply(assays_list[assay_names], t)
            assays_list <- lapply(assays_list, .makeNumpyFriendly)
        } else {
            assays_list <- rep(list(fake_mat), length(assay_names))
            names(assays_list) <- assay_names
        }
        adata$layers <- assays_list
    }

    red_dims <- as.list(reducedDims(sce))
    red_dims <- lapply(red_dims, .makeNumpyFriendly)
    adata$obsm <- red_dims

    meta_list <- metadata(sce)
    uns_list <- list()
    for (item_name in names(meta_list)) {
        item <- meta_list[[item_name]]
        tryCatch({
            # Try to convert the item using reticulate, skip if it fails
            # Capture the object output printed by reticulate
            capture.output(r_to_py(item))
            uns_list[[item_name]] <- item
        }, error = function(err) {
            warning("the '", item_name, "' item in 'metadata' cannot be ",
                    "converted to a Python type and has been skipped")
        })
    }

    adata$uns$data <- uns_list

    adata$varp <- as.list(rowPairs(sce, asSparse=TRUE))
    adata$obsp <- as.list(colPairs(sce, asSparse=TRUE))

    if (!is.null(colnames(sce))) {
        adata$obs_names <- colnames(sce)
    }

    if (!is.null(rownames(sce))) {
        adata$var_names <- rownames(sce)
    }

    adata
}

#' @importFrom DelayedArray is_sparse
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.makeNumpyFriendly <- function(x) {
    # Written originally by Charlotte Soneson in kevinrue/velociraptor.
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}
