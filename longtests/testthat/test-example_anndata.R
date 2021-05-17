library(SingleCellExperiment)

file <- system.file("extdata", "example_anndata.h5ad", package = "zellkonverter")
outfile <- tempfile(fileext = ".h5ad")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    # uns["rank_genes_groups"] doesn't convert
    expect_identical(
        names(metadata(sce)),
        c("louvain", "neighbors", "pca", "umap")
    )
    expect_identical(assayNames(sce), c("X", "counts"))
    expect_identical(
        colnames(rowData(sce)),
        c("n_counts", "highly_variable", "means", "dispersions",
          "dispersions_norm", "varm")
    )
    expect_identical(colnames(colData(sce)), "louvain")
    expect_identical(reducedDimNames(sce), c("X_pca", "X_umap"))
    expect_identical(colnames(rowData(sce)$varm), "PCs")
    expect_identical(names(colPairs(sce)), c("connectivities", "distances"))
})

test_that("Writing H5AD works", {
    sce <- readH5AD(file)
    writeH5AD(sce, outfile)
    expect_true(file.exists(outfile))
})

test_that("Round trip is as expected", {
    out <- readH5AD(outfile)

    expect_identical(dimnames(out), dimnames(sce))
    expect_identical(metadata(out), metadata(sce))
    expect_identical(assayNames(out), assayNames(sce))
    for (assay in assayNames(sce)) {
        expect_equal(assay(out, assay), assay(sce, assay))
    }
    expect_identical(reducedDims(out), reducedDims(sce))
    expect_identical(rowData(out), rowData(sce))
    expect_identical(colData(out), colData(sce))
    expect_identical(rowPairs(out), rowPairs(sce))
    expect_identical(colPairs(out), colPairs(sce))
})
