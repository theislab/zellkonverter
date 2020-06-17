#' AnnData to/from SingleCellExperiment
#'
#' Converts a Python AnnData object to or from a \linkS4class{SingleCellExperiment} object.
#'
#' @param adata A \pkg{reticulate} reference to a Python AnnData object.
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name Name of the assay to use as the primary matrix (\code{X}) of the AnnData object.
#' If `NULL`, the first assay of \code{sce} will be used by default.
#'
#' @details
#' These functions assume that an appropriate Python environment has already been loaded.
#' As such, they are largely intended for developer use, most typically inside a \pkg{basilisk} context.
#'
#' The conversion is not entirely lossless.
#' No attempt is made by \code{AnnData2SCE} to transfer the alternative Experiments from \code{sce} to an AnnData object.
#' Conversely, values in the \code{obsm} field of \code{adata} are not transferred to a SingleCellExperiment.
#'
#' @author Luke Zappia
#' 
#' @return \code{AnnData2SCE} will return a SingleCellExperiment containing the equivalent data from \code{adata}.
#'
#' \code{SCE2AnnData} will return a \pkg{reticulate} reference to an AnnData object containing the content of \code{sce}.
#'
#' @seealso
#' \code{\link{writeH5AD}} and \code{\link{readH5AD}}, for more user-friendly versions of these functions.
#' 
#' @examples
#' library(basilisk)
#' library(scRNAseq)
#' seger <- SegerstolpePancreasData()
#'
#' # If you don't know what the code below is doing,
#' # you probably shouldn't be using these functions.
#' roundtrip <- basiliskRun(env=zellkonverter:::anndata_env, fun=function(sce) {
#'      # Convert SCE to AnnData:
#'      ad <- SCE2AnnData(sce) 
#'
#'      # Maybe do some work in Python on 'ad':
#'      # BLAH BLAH BLAH
#'
#'      # Convert back to an SCE:
#'      AnnData2SCE(ad)
#' }, sce=seger)
#' 
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

#' @export
#' @rdname AnnData2SCE
SCE2AnnData <- function(sce, X_name = NULL) {

    anndata <- reticulate::import("anndata")

    if (is.null(X_name)) {
        X_name <- SummarizedExperiment::assayNames(sce)[1]
        message("using the '", X_name, "' assay as the X matrix")
    }

    X <- t(assay(sce, X_name))
    adata <- anndata$AnnData(X = .make_np_friendly(X))

    cd <- colData(sce)
    if (ncol(cd) > 0) {
        obs <- do.call(data.frame, c(as.list(cd), check.names=FALSE, stringsAsFactors=FALSE))
        adata$obs <- obs
    }

    rd <- rowData(sce)
    if (ncol(rd) > 0) {
        var <- do.call(data.frame, c(as.list(rd), check.names=FALSE, stringsAsFactors=FALSE))
        adata$var <- var
    }

    assay_names <- assayNames(sce)
    assay_names <- assay_names[!assay_names == X_name]
    if (length(assay_names) > 0) {
        assays <- assays(sce, withDimnames = FALSE)
        assays <- lapply(assays[assay_names], t)
        assays <- lapply(assays, .make_np_friendly)
        adata$layers <- assays
    }

    red_dims <- as.list(reducedDims(sce))
    red_dims <- lapply(red_dims, .make_np_friendly)
    adata$obsm <- red_dims

    adata$obs_names <- colnames(sce)
    adata$var_names <- rownames(sce)

    adata
}

#' @importFrom DelayedArray is_sparse
#' @importClassesFrom Matrix dgCMatrix
.make_np_friendly <- function(x) {
    # Written originally by Charlotte Soneson
    # in kevinrue/velociraptor.
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}
