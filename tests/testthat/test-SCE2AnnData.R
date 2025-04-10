test_that(".makeNumpyFriendly() works correctly", {
    mat <- matrix(1:50, nrow = 10, ncol = 5)

    friendly_mat <- .makeNumpyFriendly(mat, transpose = TRUE)
    expect_identical(friendly_mat, t(mat))
    expect_identical(dim(friendly_mat), rev(dim(mat)))

    friendly_mat <- .makeNumpyFriendly(mat, transpose = FALSE)
    expect_identical(friendly_mat, mat)
    expect_identical(dim(friendly_mat), dim(mat))

    sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
    friendly_sparse_mat <- .makeNumpyFriendly(sparse_mat, transpose = TRUE)
    expect_s4_class(friendly_sparse_mat, "dgRMatrix")
    expect_identical(dim(friendly_sparse_mat), rev(dim(sparse_mat)))

    friendly_sparse_mat <- .makeNumpyFriendly(sparse_mat, transpose = FALSE)
    expect_s4_class(friendly_sparse_mat, "dgCMatrix")
    expect_identical(dim(friendly_sparse_mat), dim(sparse_mat))

    delayed_mat <- DelayedArray::DelayedArray(mat)
    friendly_delayed_mat <- .makeNumpyFriendly(delayed_mat, transpose = TRUE)
    expect_identical(friendly_delayed_mat, t(mat))
    expect_identical(dim(friendly_delayed_mat), rev(dim(mat)))

    friendly_delayed_mat <- .makeNumpyFriendly(delayed_mat, transpose = FALSE)
    expect_identical(friendly_delayed_mat, mat)
    expect_identical(dim(friendly_delayed_mat), dim(mat))

    sparse_delayed_mat <- DelayedArray::DelayedArray(sparse_mat)
    friendly_sparse_delayed_mat <- .makeNumpyFriendly(sparse_delayed_mat, transpose = TRUE)
    expect_s4_class(friendly_sparse_delayed_mat, "dgRMatrix")
    expect_identical(dim(friendly_sparse_delayed_mat), rev(dim(sparse_delayed_mat)))

    friendly_sparse_delayed_mat <- .makeNumpyFriendly(sparse_delayed_mat, transpose = FALSE)
    expect_s4_class(friendly_sparse_delayed_mat, "dgCMatrix")
    expect_identical(dim(friendly_sparse_delayed_mat), dim(sparse_delayed_mat))
})
