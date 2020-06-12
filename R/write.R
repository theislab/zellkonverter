#' Write H5AD
#'
#' Write a H5AD file from a \linkS4class{SingleCellExperiment} object.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param file String containing a path to write the new `.h5ad` file.
#'
#' @details
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies.
#' This will not be performed on any subsequent run or if any other
#' \pkg{zellkonverter} function has been run prior to this one.
#'
#' @return A `NULL` is invisibly returned.
#' @author Luke Zappia
#' @seealso
#' \code{\link{readH5AD}}, to read a SingleCellExperiment file from a H5AD file.
#'
#' \code{\link{SCE2AnnData}}, for developers to create an AnnData object from a SingleCellExperiment.
#'
#' @examples
#' # Using our old friend, the Zeisel brain dataset.
#' library(scRNAseq)
#' sce <- ZeiselBrainData()
#'
#' # Writing to a H5AD file.
#' temp <- tempfile(fileext='.h5ad')
#' writeH5AD(sce, temp)
#'
#' @export
writeH5AD <- function(sce, file) {
    proc <- basilisk::basiliskStart(anndata_env)
    on.exit(basilisk::basiliskStop(proc))

    file <- path.expand(file)

    basilisk::basiliskRun(proc, .H5ADwriter, sce=sce, file=file)

    invisible(NULL)
}

.H5ADwriter <- function(sce, file) {
    anndata <- reticulate::import("anndata")
    adata <- SCE2AnnData(sce)
    adata$write_h5ad(file)
}
