library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30578574")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X"),
    colData = c("paul15_clusters", "n_counts_all", "louvain"),
    rowData = c("n_counts", "mean", "std"),
    metadata = c("draw_graph", "iroot", "louvain", "louvain_sizes", "neighbors",
                 "paga", "pca"),
    redDim = c("X_draw_graph_fa", "X_pca"),
    varm = c("PCs"),
    colPairs = c("connectivities", "distances")
)

missing <- list(
    metadata = c("rank_genes_groups")
)

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
