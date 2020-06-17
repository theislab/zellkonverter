# This tests the readH5AD function (and by implication, SCE2AnnData).
# library(testthat); library(zellkonverter); source("test-read.R")

library(SummarizedExperiment)
file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with HDF5Arrays", {
    sce <- readH5AD(file, use.hdf5=TRUE)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_s4_class(DelayedArray::seed(assay(sce)), "HDF5ArraySeed")
})

test_that("readH5AD works in a separate process", {
    oldshare <- basilisk::getBasiliskShared()
    basilisk::setBasiliskShared(FALSE)
    oldfork <- basilisk::getBasiliskFork()
    basilisk::setBasiliskFork(FALSE)

    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    basilisk::setBasiliskShared(oldshare)
    basilisk::setBasiliskFork(oldfork)
})

