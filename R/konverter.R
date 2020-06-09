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
