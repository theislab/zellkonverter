#' Write H5AD
#'
#' Write a H5AD file from a SingleCellExperiment
#'
#' @param sce SingleCellExperiment object
#' @param file Path to write new `.h5ad`
#'
#' @return SingleCellExperiment
#' @export
writeH5AD <- function(sce, file) {
    proc <- basilisk::basiliskStart(anndata_env)
    on.exit(basilisk::basiliskStop(proc))

    file <- path.expand(file)

    basilisk::basiliskRun(proc, .H5ADwriter, sce=sce, file=file)
}

.H5ADwriter <- function(sce, file) {
    anndata <- reticulate::import("anndata")
    adata <- SCE2AnnData(sce)
    adata$write_h5ad(file)
}
