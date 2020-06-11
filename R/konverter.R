#' AnnData to SingleCellExperiment
#'
#' Converts a Python AnnData object to a SingleCellExperiment object.
#'
#' @param adata Reference to a Python AnnData object
#'
#' @return SingleCellExperiment
#' @export
AnnData2SCE <- function(adata) {

    py_builtins <- reticulate::import_builtins()

    x_mat <- t(adata$X)
    colnames(x_mat) <- adata$obs_names$to_list()
    rownames(x_mat) <- adata$var_names$to_list()

    assays_list <- list(X = x_mat)
    layer_names <- names(py_builtins$dict(adata$layers))

    if (length(layer_names) > 0) {
        for (layer_name in layer_names) {
            assays_list[[layer_name]] <- t(adata$layers$get(layer_name))
        }
    }

    SingleCellExperiment::SingleCellExperiment(
        assays      = assays_list,
        rowData     = adata$var,
        colData     = adata$obs,
        reducedDims = py_builtins$dict(adata$obsm)
    )
}

#' SingleCellExperiment to AnnData
#'
#' Converts a SingleCellExperiment object a Python AnnData object.
#'
#' @param sce SingleCellExperiment object
#' @param X_name Name of the assay to use as the X matrix of the AnnData object.
#' If `NULL` the first assay will be used.
#'
#' @return Reference to a Python AnnData object
#' @export
SCE2AnnData <- function(sce, X_name = NULL) {

    anndata <- reticulate::import("anndata")

    if (is.null(X_name)) {
        X_name <- SummarizedExperiment::assayNames(sce)[1]
        message("Using the '", X_name, "' assay as the X matrix")
    }

    X <- t(SummarizedExperiment::assay(sce, X_name))

    obs <- as.data.frame(SummarizedExperiment::colData(sce))
    var <- as.data.frame(SummarizedExperiment::rowData(sce))

    adata <- anndata$AnnData(X = X)

    if (ncol(obs) > 0) {
        adata$obs <- obs
    }

    if (ncol(var) > 0) {
        adata$var <- var
    }

    assay_names <- SummarizedExperiment::assayNames(sce)
    assay_names <- assay_names[!assay_names == X_name]

    if (length(assay_names) > 0) {
        assays <- SummarizedExperiment::assays(sce, withDimnames = FALSE)
        adata$layers <- lapply(assays[assay_names], t)
    }

    adata$obsm <- as.list(SingleCellExperiment::reducedDims(sce))

    adata$obs_names <- colnames(sce)
    adata$var_names <- rownames(sce)

    adata
}
