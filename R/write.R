#' Write H5AD
#'
#' Write a H5AD file from a \linkS4class{SingleCellExperiment} object.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param file String containing a path to write the new `.h5ad` file.
#' @param skip_assays Logical scalar indicating whether assay matrices should
#' be ignored when writing to `file`.
#'
#' @details
#' Setting `skip_assays=TRUE` can occasionally be useful if the matrices in
#' `sce` are stored in a format that is not amenable for efficient conversion
#' to a **numpy**-compatible format. In such cases, it can be better to create
#' an empty placeholder dataset in `file` and fill it in R afterwards.
#'
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
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(scRNAseq)
#'     sce <- ZeiselBrainData()
#'
#'     # Writing to a H5AD file
#'     temp <- tempfile(fileext = '.h5ad')
#'     writeH5AD(sce, temp)
#' }
#'
#' @export
#' @importFrom basilisk basiliskRun
writeH5AD <- function(sce, file, skip_assays = FALSE) {
    file <- path.expand(file)
    basiliskRun(env = anndata_env, fun = .H5ADwriter, sce = sce, file = file, skip_assays = skip_assays)
    invisible(NULL)
}

#' @importFrom reticulate import
.H5ADwriter <- function(sce, file, skip_assays) {
    anndata <- import("anndata")
    adata <- SCE2AnnData(sce, skip_assays = skip_assays)
    adata$write_h5ad(file)
}
