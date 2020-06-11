#' Read H5AD
#'
#' Reads a H5AD file and returns a \linkS4class{SingleCellExperiment} object.
#'
#' @param file String containing a path to a `.h5ad` file.
#'
#' @details
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies.
#' This will not be performed on any subsequent run or if any other
#' \pkg{zellkonverter} function has been run prior to this one.
#'
#' @return A \linkS4class{SingleCellExperiment} object is returned.
#'
#' @examples
#' # Re-using the example from writeH5AD.
#' example(writeH5AD, echo=FALSE)
#'
#' # Reading into a SingleCellExperiment.
#' sce2 <- readH5AD(temp)
#'
#' @export
readH5AD <- function(file) {
    proc <- basilisk::basiliskStart(anndata_env)
    on.exit(basilisk::basiliskStop(proc))

    file <- path.expand(file)

    basilisk::basiliskRun(proc, .H5ADreader, file=file)
}

.H5ADreader <- function(file) {
    anndata <- reticulate::import("anndata")
    adata <- anndata$read_h5ad(file)
    AnnData2SCE(adata)
}
