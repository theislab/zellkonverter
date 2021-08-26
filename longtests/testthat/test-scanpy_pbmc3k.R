library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://ndownloader.figshare.com/files/30462915")
outfile <- tempfile(fileext = ".h5ad")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), c("X"))

    expect_identical(
        colnames(colData(sce)),
        c("n_genes", "n_genes_by_counts", "total_counts", "total_counts_mt",
          "pct_counts_mt", "leiden")
    )

    expect_identical(
        colnames(rowData(sce)),
        c("gene_ids", "n_cells", "mt", "n_cells_by_counts", "mean_counts",
          "pct_dropout_by_counts", "total_counts", "highly_variable", "means",
          "dispersions", "dispersions_norm", "mean", "std", "varm")
    )

    expect_identical(
        names(metadata(sce)),
        c("hvg", "leiden", "neighbors", "pca", "umap")
    )

    # Known issue with converting rank_genes_groups.
    # This should fail once that has been fixed.
    expect_true(
        !("rank_genes_groups" %in% names(metadata(sce)))
    )

    expect_identical(
        reducedDimNames(sce),
        c("X_pca", "X_umap")
    )

    expect_identical(
        colnames(rowData(sce)$varm),
        c("PCs")
    )

    expect_identical(
        names(colPairs(sce)),
        c("connectivities", "distances")
    )
})

test_that("Writing H5AD works", {
    sce <- readH5AD(file)
    writeH5AD(sce, outfile)
    expect_true(file.exists(outfile))
})

test_that("Round trip is as expected", {
    out <- readH5AD(outfile)

    expect_identical(dimnames(out), dimnames(sce))
    if (length(metadata(sce)) > 0) {
        expect_identical(metadata(out), metadata(sce))
    }
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
