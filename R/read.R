#' Read H5AD
#'
#' Reads a H5AD file and returns a \linkS4class{SingleCellExperiment} object.
#'
#' @param file String containing a path to a `.h5ad` file.
#' @param use_hdf5 Logical scalar indicating whether assays should be
#' loaded as HDF5-based matrices from the **HDF5Array** package.
#'
#' @details
#' Setting `use_hdf5 = TRUE` allows for very large datasets to be efficiently
#' represented on machines with little memory. However, this comes at the cost
#' of access speed as data needs to be fetched from the HDF5 file upon request.
#'
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies. This will not be performed on
#' any subsequent run or if any other **zellkonverter** function has been run
#' prior to this one.
#'
#' @return A \linkS4class{SingleCellExperiment} object is returned.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
#' sce <- readH5AD(file)
#' class(assay(sce))
#'
#' sce2 <- readH5AD(file, use_hdf5 = TRUE)
#' class(assay(sce2))
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @seealso
#' [`writeH5AD()`], to write a \linkS4class{SingleCellExperiment} object to a
#' H5AD file.
#'
#' [`AnnData2SCE()`], for developers to convert existing AnnData instances to a
#' \linkS4class{SingleCellExperiment}.
#'
#' @export
#' @importFrom basilisk basiliskRun 
readH5AD <- function(file, use_hdf5 = FALSE) {
    file <- path.expand(file)

    basiliskRun(
        env    = anndata_env,
        fun    = .H5ADreader,
        file   = file,
        backed = use_hdf5
    )
}

#' @importFrom reticulate import
.H5ADreader <- function(file, backed = FALSE) {
    anndata <- import("anndata")
    adata <- anndata$read_h5ad(file, backed = if (backed) "r" else FALSE)
    AnnData2SCE(adata, hdf5_backed = backed)
}
