# This tests the writeH5AD function (and by implication, AnnData2SCE).
# library(testthat); library(zellkonverter); source("test-write.R")

library(scRNAseq)
sce <- ZeiselBrainData()
reducedDim(sce, "WHEE") <- matrix(runif(ncol(sce) * 10), ncol = 10)

test_that("writeH5AD works as expected", {
    temp <- tempfile(fileext = '.h5ad')
    writeH5AD(sce, temp)

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
