#' Read H5AD
#'
#' Reads a H5AD file and returns a \linkS4class{SingleCellExperiment} object.
#'
#' @param file String containing a path to a `.h5ad` file.
#' @param use.hdf5 Logical scalar indicating whether assays should be
#' loaded as HDF5-based matrices from the \pkg{HDF5Array} package.
#'
#' @details
#' Setting \code{use.hdf5=TRUE} allows for very large datasets to be
#' efficiently represented on machines with little memory.
#' However, this comes at the cost of access speed as data needs to 
#' be fetched from the HDF5 file upon request. 
#'
#' When first run, this function will instantiate a conda environment
#' containing all of the necessary dependencies.
#' This will not be performed on any subsequent run or if any other
#' \pkg{zellkonverter} function has been run prior to this one.
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
#' sce2 <- readH5AD(file, use.hdf5=TRUE)
#' class(assay(sce2))
#'
#' @author Luke Zappia
#' @seealso
#' \code{\link{writeH5AD}}, to write a SingleCellExperiment file to a H5AD file.
#'
#' \code{\link{AnnData2SCE}}, for developers to convert existing AnnData instances to a SingleCellExperiment.
#'
#' @export
readH5AD <- function(file, use.hdf5=FALSE) {
    file <- path.expand(file)

    # We set shared=!use.hdf5 because AnnData opens a blocking r+ connection to
    # the HDF5 file that is really hard to shut down via reticulate.
    output <- basilisk::basiliskRun(env=anndata_env, shared=!use.hdf5,
        fun=.H5ADreader, file=file, backed=use.hdf5)

    if (use.hdf5) {
        SummarizedExperiment::assay(output, "X", withDimnames=FALSE) <- HDF5Array::HDF5Array(file, "X")
        for (i in setdiff(SummarizedExperiment::assayNames(output), "X")) {
            SummarizedExperiment::assay(output, i, withDimnames=FALSE) <- HDF5Array::HDF5Array(file, file.path("layers", i))
        }
    }

    output
}

.H5ADreader <- function(file, backed=FALSE) {
    anndata <- reticulate::import("anndata")
    adata <- anndata$read_h5ad(file, backed=backed)
    output <- AnnData2SCE(adata, skip.assays=backed)
}
