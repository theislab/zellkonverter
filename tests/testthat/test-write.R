# This tests the writeH5AD function (and by implication, AnnData2SCE).
# library(testthat); library(zellkonverter); source("test-write.R")

library(scRNAseq)
sce <- ZeiselBrainData()
reducedDim(sce, "WHEE") <- matrix(runif(ncol(sce) * 10), ncol = 10)

test_that("writeH5AD works as expected", {
    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp)

    expect_identical(dimnames(out), dimnames(sce))
    expect_equal(assay(out), assay(sce))
    expect_identical(reducedDims(out), reducedDims(sce))

    # Need to coerce the factors back to strings.
    row_data <- rowData(out)
    for (i in seq_len(ncol(row_data))) {
        if (is.factor(row_data[[i]])) {
            row_data[[i]] <- as.character(row_data[[i]])
        }
    }
    expect_identical(row_data, rowData(sce))

    col_data <- colData(out)
    for (i in seq_len(ncol(col_data))) {
        if (is.factor(col_data[[i]])) {
            col_data[[i]] <- as.character(col_data[[i]])
        }
    }
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with sparse matrices", {
    mat <- assay(sce)
    counts(sce) <- as(mat, "dgCMatrix")
    logcounts(sce) <- counts(sce) * 10
    assay(sce, "random") <- mat # throwing in a dense matrix in a mixture.

    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp)

    expect_identical(counts(sce), assay(out, "X"))
    expect_identical(logcounts(sce), logcounts(out))
    # expect_identical() was failing on Windows for some reason...
    expect_equal(assay(sce, "random"), assay(out, "random"))
})

test_that("writeH5AD works with assay skipping", {
    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp, skip_assays = TRUE)
    expect_true(file.exists(temp))

    out <- HDF5Array::HDF5Array(temp, "X/data")
    expect_identical(sum(out), 0) # it's empty!
})

test_that("writeH5AD works with X_name", {
    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp, X_name = "counts")
    expect_true(file.exists(temp))

    out <- readH5AD(temp)
    expect_equal(assay(out, "X"), assay(sce, "counts"))
})

test_that("writeH5AD works in a separate process", {
    oldshare <- basilisk::getBasiliskShared()
    basilisk::setBasiliskShared(FALSE)
    oldfork <- basilisk::getBasiliskFork()
    basilisk::setBasiliskFork(FALSE)

    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))

    basilisk::setBasiliskShared(oldshare)
    basilisk::setBasiliskFork(oldfork)
})
