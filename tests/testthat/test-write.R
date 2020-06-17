# This tests the writeH5AD function (and by implication, AnnData2SCE).
# library(testthat); library(zellkonverter); source("test-write.R")

library(scRNAseq)
sce <- ZeiselBrainData()
reducedDim(sce, "WHEE") <- matrix(runif(ncol(sce)*10), ncol=10)

test_that("writeH5AD works as expected", {
    temp <- tempfile(fileext='.h5ad')
    writeH5AD(sce, temp)

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp)

    expect_identical(dimnames(out), dimnames(sce))
    expect_equal(assay(out), assay(sce))
    expect_identical(reducedDims(out), reducedDims(sce))

    # Need to coerce the factors back to strings. 
    rd <- rowData(out)
    for (i in seq_len(ncol(rd))) {
        if (is.factor(rd[[i]])) {
            rd[[i]] <- as.character(rd[[i]])
        }
    }
    expect_identical(rd, rowData(sce))

    cd <- colData(out)
    for (i in seq_len(ncol(cd))) {
        if (is.factor(cd[[i]])) {
            cd[[i]] <- as.character(cd[[i]])
        }
    }
    expect_identical(cd, colData(sce))
})

test_that("writeH5AD works in a separate process", {
    oldshare <- basilisk::getBasiliskShared()
    basilisk::setBasiliskShared(FALSE) 
    oldfork <- basilisk::getBasiliskFork()
    basilisk::setBasiliskFork(FALSE)

    temp <- tempfile(fileext='.h5ad')
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))

    basilisk::setBasiliskShared(oldshare)
    basilisk::setBasiliskFork(oldfork)
})
