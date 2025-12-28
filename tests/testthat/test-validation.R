file <- system.file("extdata", "example_anndata.h5ad",
    package = "zellkonverter"
)
sce <- readH5AD(file)

names <- list(
    assays = c("X", "counts"),
    colData = "louvain",
    rowData = c(
        "n_counts", "highly_variable", "means", "dispersions",
        "dispersions_norm"
    ),
    metadata = c("louvain", "neighbors", "pca", "rank_genes_groups", "umap"),
    redDim = c("X_pca", "X_umap"),
    varm = "PCs",
    colPairs = c("connectivities", "distances")
)

missing <- list()

test_that("validateH5ADSCE works", {
    validateH5ADSCE(sce, names, missing)
    expect_error(
        validateH5ADSCE(sce, names, list(varm = "PCs")),
        "Expected varm names missing to be TRUE"
    )
})

test_that("expectSCE works", {
    expectSCE(sce, sce)
})
