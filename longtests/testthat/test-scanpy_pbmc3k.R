library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://ndownloader.figshare.com/files/30462915")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X"),
    colData = c("n_genes", "n_genes_by_counts", "total_counts",
                "total_counts_mt", "pct_counts_mt", "leiden"),
    rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts", "mean_counts",
                "pct_dropout_by_counts", "total_counts", "highly_variable",
                "means", "dispersions", "dispersions_norm", "mean", "std"),
    metadata = c("hvg", "leiden", "neighbors", "pca", "rank_genes_groups",
                 "umap"),
    redDim = c("X_pca", "X_umap"),
    varm = c("PCs"),
    colPairs = c("connectivities", "distances")
)

missing <- list()

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
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
    out <- readH5AD(outfile)
    expectSCE(out, sce)
})
