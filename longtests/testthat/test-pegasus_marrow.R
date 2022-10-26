library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30682400")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X"),
    colData = c("n_genes", "Channel", "n_counts", "percent_mito", "scale",
                "Group", "louvain_labels", "anno"),
    rowData = c("featureid", "n_cells", "percent_cells", "robust",
                "highly_variable_features", "mean", "var", "hvf_loess",
                "hvf_rank"),
    metadata = c("Channels", "Groups", "PCs", "W_diffmap", "W_pca_harmony",
                 "c2gid", "diffmap_evals", "diffmap_knn_distances",
                 "diffmap_knn_indices", "genome", "gncells",
                 "louvain_resolution", "modality", "ncells", "norm_count",
                 "pca", "pca_features", "pca_harmony_knn_distances",
                 "pca_harmony_knn_indices", "stdzn_max_value", "stdzn_mean",
                 "stdzn_std"),
    redDim = c("X_diffmap", "X_fle", "X_pca", "X_pca_harmony", "X_phi",
               "X_tsne", "X_umap"),
    varm = c("de_res", "gmeans", "gstds", "means", "partial_sum")
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
