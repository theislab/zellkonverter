library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://ndownloader.figshare.com/files/24539828")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), c("X", "counts"))

    expect_identical(
        colnames(colData(sce)),
        c("tech", "celltype", "size_factors")
    )
})

test_that("Writing H5AD works", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))
})

test_that("Round trip is as expected", {
    out <- readH5AD(temp)

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
