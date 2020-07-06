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
    sce <- readH5AD(file, use_hdf5 = TRUE)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_s4_class(DelayedArray::seed(assay(sce)), "HDF5ArraySeed")

    ref <- readH5AD(file)
    expect_identical(as.matrix(assay(ref)), as.matrix(assay(sce)))

    # Properly sleeps to wait for the process to shut down. 
    expect_s4_class(sce <- readH5AD(file, use_hdf5 = TRUE), "SingleCellExperiment")
})

test_that("Reading H5AD works with a mixture of sparse and HDF5Arrays", {
    sce <- readH5AD(file)
    assay(sce, "more") <- as(assay(sce, "X"), "dgCMatrix")

    temp <- tempfile(fileext=".h5ad")
    writeH5AD(sce, temp)

    backed <- readH5AD(temp, use_hdf5 = TRUE)
    expect_s4_class(DelayedArray::seed(assay(backed)), "HDF5ArraySeed")
    expect_s4_class(assay(backed, "more"), "dgCMatrix")
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
