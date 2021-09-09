library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30639279")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X"),
    colData = c("in_tissue", "array_row", "array_col", "n_genes_by_counts",
                "log1p_n_genes_by_counts", "total_counts", "log1p_total_counts",
                "pct_counts_in_top_50_genes", "pct_counts_in_top_100_genes",
                "pct_counts_in_top_200_genes", "pct_counts_in_top_500_genes",
                "total_counts_MT", "log1p_total_counts_MT", "pct_counts_MT",
                "n_counts", "leiden", "cluster", "features_summary_cluster",
                "features_histogram_cluster", "features_texture_cluster"),
    rowData = c("gene_ids", "feature_types", "genome", "MT",
                "n_cells_by_counts", "mean_counts", "log1p_mean_counts",
                "pct_dropout_by_counts", "total_counts", "log1p_total_counts",
                "n_cells", "highly_variable", "highly_variable_rank", "means",
                "variances", "variances_norm"),
    metadata = c("cluster_co_occurrence", "cluster_colors", "cluster_ligrec",
                 "cluster_nhood_enrichment", "hvg", "leiden", "leiden_colors",
                 "moranI", "neighbors", "pca", "spatial", "spatial_neighbors",
                 "umap"),
    redDim = c("X_pca", "X_umap", "features", "features_context",
               "features_lowres", "features_orig", "features_segmentation",
               "spatial"),
    varm = c("PCs"),
    colPairs = c("connectivities", "distances", "spatial_connectivities",
                 "spatial_distances")
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
