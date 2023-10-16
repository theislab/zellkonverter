# This tests the writeH5AD function (and by implication, AnnData2SCE).
library(scRNAseq)

sce <- ZeiselBrainData()
reducedDim(sce, "WHEE") <- matrix(runif(ncol(sce) * 10), ncol = 10)

test_that("writeH5AD works as expected", {
    temp <- tempfile(fileext = ".h5ad")
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
    names(col_data) <- names(colData(sce))
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with version 0.9.2", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, version = "0.9.2")
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp, version = "0.9.2")

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
    names(col_data) <- names(colData(sce))
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with version 0.8.0", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, version = "0.8.0")
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp, version = "0.8.0")

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
    names(col_data) <- names(colData(sce))
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with version 0.7.6", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, version = "0.7.6")
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp, version = "0.7.6")

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
    names(col_data) <- names(colData(sce))
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with verbose=TRUE", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, verbose = TRUE)
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
    names(col_data) <- names(colData(sce))
    expect_identical(col_data, colData(sce))
})

test_that("writeH5AD works as expected with sparse matrices", {
    sparse_sce <- sce
    mat <- assay(sparse_sce)
    counts(sparse_sce) <- as(mat, "CsparseMatrix")
    logcounts(sparse_sce) <- counts(sparse_sce) * 10
    assay(sparse_sce, "random") <- mat # throwing in a dense matrix in a mixture.

    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sparse_sce, temp)
    expect_true(file.exists(temp))

    # Reading it back out again. Hopefully we didn't lose anything important.
    out <- readH5AD(temp, X_name = "X")

    expect_identical(counts(sparse_sce), assay(out, "X"))
    expect_identical(logcounts(sparse_sce), logcounts(out))
    # expect_identical() was failing on Windows for some reason...
    expect_equal(assay(sparse_sce, "random"), assay(out, "random"))
})

test_that("writeH5AD works with assay skipping", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, skip_assays = TRUE)
    expect_true(file.exists(temp))

    out <- HDF5Array::HDF5Array(temp, "X/data")
    expect_identical(sum(out), 0) # it's empty!
})

test_that("writeH5AD works with X_name", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, X_name = "counts")
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")
    expect_equal(assay(out, "X"), assay(sce, "counts"))
})

test_that("writeH5AD works in a separate process", {
    oldshare <- basilisk::getBasiliskShared()
    basilisk::setBasiliskShared(FALSE)
    oldfork <- basilisk::getBasiliskFork()
    basilisk::setBasiliskFork(FALSE)

    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp)
    expect_true(file.exists(temp))

    basilisk::setBasiliskShared(oldshare)
    basilisk::setBasiliskFork(oldfork)
})

test_that("writeH5AD DelayedArray X works", {
    delayed_sce <- sce
    counts(delayed_sce) <- DelayedArray::DelayedArray(counts(delayed_sce))

    temp <- tempfile(fileext = ".h5ad")

    writeH5AD(delayed_sce, temp, X_name = "counts")
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")

    # Identical fail on Windows for some reason
    expect_equal(counts(sce), assay(out, "X"))
})

test_that("writeH5AD sparse DelayedArray X works", {
    delayed_sce <- sce
    sparse_counts <- as(counts(delayed_sce), "CsparseMatrix")
    counts(delayed_sce) <- DelayedArray::DelayedArray(sparse_counts)

    temp <- tempfile(fileext = ".h5ad")

    writeH5AD(delayed_sce, temp, X_name = "counts")
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")

    # Sparse DelayedArrays are currently coerced into memory
    # This expectation will need to be changed once that is fixed
    expect_identical(sparse_counts, assay(out, "X"))
})

test_that("writeH5AD DelayedArray layer works", {
    delayed_sce <- sce
    assay(delayed_sce, "layer") <- DelayedArray::DelayedArray(
        counts(delayed_sce)
    )

    temp <- tempfile(fileext = ".h5ad")

    writeH5AD(delayed_sce, temp)
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")

    # Identical fails on Windows for some reason
    expect_equal(counts(sce), assay(out, "layer"))
})

test_that("writeH5AD works with colData list columns", {
    list_sce <- sce
    colData(list_sce)$ListCol <- lapply(seq_len(ncol(list_sce)), function(x) {
        sample(LETTERS, 2)
    })

    temp <- tempfile(fileext = ".h5ad")

    expect_warning(writeH5AD(list_sce, temp), "columns are not atomic")
    expect_true(file.exists(temp))

    # Knowing what comes back is hard so just check there is something
    out <- readH5AD(temp, X_name = "X")
    expect_true("ListCol" %in% names(metadata(out)$.colData))
})

test_that("writeH5AD works with rowData list columns", {
    list_sce <- sce
    rowData(list_sce)$ListCol <- lapply(seq_len(nrow(list_sce)), function(x) {
        sample(LETTERS, 2)
    })

    temp <- tempfile(fileext = ".h5ad")

    expect_warning(writeH5AD(list_sce, temp), "columns are not atomic")
    expect_true(file.exists(temp))

    # Knowing what comes back is hard so just check there is something
    out <- readH5AD(temp, X_name = "X")
    expect_true("ListCol" %in% names(metadata(out)$.rowData))
})

test_that("writeH5AD works with gzip compression", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, X_name = "counts", compression = "gzip")
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")
    expect_equal(assay(out, "X"), assay(sce, "counts"))
})

test_that("writeH5AD works with lzf compression", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, X_name = "counts", compression = "lzf")
    expect_true(file.exists(temp))

    out <- readH5AD(temp, X_name = "X")
    expect_equal(assay(out, "X"), assay(sce, "counts"))
})

test_that("Skipping slot conversion works", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, assays = FALSE, colData = FALSE, rowData = FALSE,
              varm = FALSE, reducedDims = FALSE, metadata = FALSE,
              colPairs = FALSE, rowPairs = FALSE)

    out <- readH5AD(temp, X_name = "X")

    expect_identical(assayNames(out), "X")
    expect_identical(metadata(out), list(X_name = "counts"))
    expect_equal(ncol(rowData(out)), 0)
    expect_equal(ncol(colData(out)), 0)
    expect_equal(length(reducedDims(out)), 0)
    expect_equal(length(rowPairs(out)), 0)
    expect_equal(length(colPairs(out)), 0)
})

test_that("Selective DF conversion works", {
    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, assays = FALSE, colData = "tissue")

    out <- readH5AD(temp, X_name = "X")

    expect_identical(names(colData(out)), "tissue")
})
