#' Write H5AD
#'
#' Write a H5AD file from a \linkS4class{SingleCellExperiment} object.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param file String containing a path to write the new `.h5ad` file.
#'
#' @details
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies. This will not be performed on
#' any subsequent run or if any other
#' **zellkonverter** function has been run prior to this one.
#'
#' @return A `NULL` is invisibly returned.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @seealso
#' [`readH5AD()`], to read a \linkS4class{SingleCellExperiment} file from a H5AD
#' file.
#'
#' [`SCE2AnnData()`], for developers to create an AnnData object from a
#' \linkS4class{SingleCellExperiment}.
#'
#' @examples
#' # Using the Zeisel brain dataset
#' library(scRNAseq)
#' sce <- ZeiselBrainData()
#'
#' # Writing to a H5AD file
#' temp <- tempfile(fileext = '.h5ad')
#' writeH5AD(sce, temp)
#'
#' @export
#' @importFrom basilisk basiliskRun
writeH5AD <- function(sce, file) {
    file <- path.expand(file)
    basiliskRun(env = anndata_env, fun = .H5ADwriter, sce = sce, file = file)
    invisible(NULL)
}

#' @importFrom reticulate import
.H5ADwriter <- function(sce, file) {
    anndata <- import("anndata")
    adata <- SCE2AnnData(sce)
    adata$write_h5ad(file)
}
