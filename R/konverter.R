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
#' For `AnnData2SCE()`, an error is raised if there is no corresponding R format
#' for a matrix in the AnnData object. If `skip_assays = TRUE`, no error is
#' given and empty sparse matrices are created for each assay. The user is
#' expected to fill in the assays on the R side, see [`readH5AD()`] for an
#' example.
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
#' any assays in `sce`, replacing them with empty sparse matrices instead.
#'
#' @export
#' @importFrom Matrix t sparseMatrix
#' @importFrom methods selectMethod is
AnnData2SCE <- function(adata, skip_assays = FALSE) {
    py_builtins <- reticulate::import_builtins()

    if (!skip_assays) {
        x_mat <- adata$X
        t_FUN <- selectMethod("t", signature = class(x_mat), optional = TRUE)
        if (is.null(t_FUN)) {
            stop("assay matrices do not support transposition")
        }
        x_mat <- t_FUN(x_mat)
    } else {
        dims <- unlist(adata$X$shape)
        fake_mat <- sparseMatrix(
            i    = integer(0),
            j    = integer(0),
            x    = numeric(0),
            dims = rev(dims)
        )
        x_mat <- fake_mat
    }

    colnames(x_mat) <- adata$obs_names$to_list()
    rownames(x_mat) <- adata$var_names$to_list()

    assays_list <- list(X = x_mat)
    layer_names <- names(py_builtins$dict(adata$layers))

    if (length(layer_names) > 0) {
        for (layer_name in layer_names) {
            if (!skip_assays) {
                layer_mat <- t(adata$layers$get(layer_name))
            } else {
                layer_mat <- fake_mat
            }
            assays_list[[layer_name]] <- layer_mat
        }
    }

    uns_data <- adata$uns$data

    meta_list <- list()
    for (item_name in names(uns_data)) {
        item <- uns_data[[item_name]]
        if (!(is(item, "python.builtin.object"))) {
            meta_list[[item_name]] <- item
        } else {
            warning("the '", item_name, "' item in 'uns' cannot be converted ",
                    "to an R object and has been skipped")
        }
    }

    varp_list <- lapply(py_builtins$dict(adata$varp), function(v) {v$todense()})
    obsp_list <- lapply(py_builtins$dict(adata$obsp), function(v) {v$todense()})

    row_data <- S4Vectors::DataFrame(adata$var)
    varm_list <- py_builtins$dict(adata$varm)
    if (length(varm_list) > 0) {
        # Create an empty DataFrame with the correct number of rows
        varm_df <- S4Vectors::DataFrame(matrix(, nrow = adata$n_vars, ncol = 0))
        for (varm_name in names(varm_list)) {
            varm_df[[varm_name]] <- varm_list[[varm_name]]
        }
        row_data$varm <- varm_df
    }

    SingleCellExperiment::SingleCellExperiment(
        assays      = assays_list,
        rowData     = row_data,
        colData     = adata$obs,
        reducedDims = py_builtins$dict(adata$obsm),
        metadata    = meta_list,
        rowPairs    = varp_list,
        colPairs    = obsp_list
    )
}

#' @rdname AnnData-Conversion
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name Name of the assay to use as the primary matrix (`X`) of the
#' AnnData object. If `NULL`, the first assay of `sce` will be used by default.
#'
#' @export
#' @importFrom utils capture.output
SCE2AnnData <- function(sce, X_name = NULL) {

    anndata <- reticulate::import("anndata")

    if (is.null(X_name)) {
        if (length(assays(sce)) == 0) {
            stop("'sce' does not contain any assays")
        }
        X_name <- SummarizedExperiment::assayNames(sce)[1]
        message("Note: using the '", X_name, "' assay as the X matrix")
    }

    X <- t(assay(sce, X_name))
    adata <- anndata$AnnData(X = .makeNumpyFriendly(X))

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
        assays_list <- assays(sce, withDimnames = FALSE)
        assays_list <- lapply(assays_list[assay_names], t)
        assays_list <- lapply(assays_list, .makeNumpyFriendly)
        adata$layers <- assays_list
    }

    red_dims <- as.list(reducedDims(sce))
    red_dims <- lapply(red_dims, .makeNumpyFriendly)
    adata$obsm <- red_dims

    meta_list <- S4Vectors::metadata(sce)
    uns_list <- list()
    for (item_name in names(meta_list)) {
        item <- meta_list[[item_name]]
        tryCatch({
            # Try to convert the item using reticulate, skip if it fails
            # Capture the object output printed by reticulate
            capture.output(reticulate::r_to_py(item))
            uns_list[[item_name]] <- item
        }, error = function(err) {
            warning("the '", item_name, "' item in 'metadata' cannot be ",
                    "converted to a Python type and has been skipped")
        })
    }

    adata$uns$data <- uns_list

    adata$varp <- as.list(SingleCellExperiment::rowPairs(sce, asSparse=TRUE))
    adata$obsp <- as.list(SingleCellExperiment::colPairs(sce, asSparse=TRUE))

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
