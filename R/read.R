#' Read H5AD
#'
#' Reads a H5AD file and returns a SingleCellExperiment
#'
#' @param file Path to a `.h5ad` file
#'
#' @return SingleCellExperiment
#' @export
readH5AD <- function(file) {
    proc <- basilisk::basiliskStart(anndata_env)
    on.exit(basilisk::basiliskStop(proc))

    file <- path.expand(file)

    basilisk::basiliskRun(proc, function(file) {
        anndata <- reticulate::import("anndata")
        adata <- anndata$read_h5ad(file)
        sce <- AnnData2SCE(adata)
    }, file=file)
}
