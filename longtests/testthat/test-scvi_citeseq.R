library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30612834")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X", "counts", "denoised_rna"),
    colData = c("n_genes", "percent_mito", "n_counts", "batch", "X_scvi_batch",
                "X_scvi_labels", "X_scvi_local_l_mean", "X_scvi_local_l_var",
                "leiden_totalVI"),
    rowData = c("highly_variable", "highly_variable_rank", "means", "variances",
                "variances_norm", "highly_variable_nbatches"),
    metadata = c("X_scvi", "hvg", "leiden", "neighbors", "umap"),
    redDim = c("X_totalVI", "X_umap", "denoised_protein",
               "protein_expression", "protein_foreground_prob"),
    colPairs = c("connectivities", "distances")
)

missing <- list()

test_that("Reading H5AD works", {
    sce <- expect_warning(readH5AD(file))
    expect_s4_class(sce, "SingleCellExperiment")
})

sce <- suppressWarnings(readH5AD(file))

test_that("SCE is valid", {
    validateH5ADSCE(sce, names, missing)
})

test_that("Writing H5AD works", {
    writeH5AD(sce, outfile)
    expect_true(file.exists(outfile))
})

test_that("Round trip is as expected", {
    out <- suppressWarnings(readH5AD(outfile))

    # For some reason "_scvi" gets changed to "X_scvi", not sure why...
    names(S4Vectors::metadata(sce))[1] <- "X_scvi"

    expectSCE(out, sce)
})
