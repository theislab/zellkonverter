# This tests the readH5AD function (and by implication, SCE2AnnData).
library(SummarizedExperiment)
file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
file_v08 <- system.file("extdata", "krumsiek11_augmented_v0-8.h5ad", package = "zellkonverter")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with version 0.9.2", {
    sce <- readH5AD(file, version = "0.9.2")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with version 0.8.0", {
    sce <- readH5AD(file, version = "0.8.0")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with version 0.7.6", {
    sce <- readH5AD(file, version = "0.7.6")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with verbose=TRUE", {
    sce <- readH5AD(file, verbose = TRUE)
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
    expect_s4_class(
        sce <- readH5AD(file, use_hdf5 = TRUE),
        "SingleCellExperiment"
    )
})

test_that("Reading H5AD works with a mixture of sparse and HDF5Arrays", {
    sce <- readH5AD(file)
    assay(sce, "more") <- as(assay(sce, "X"), "CsparseMatrix")

    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp)

    backed <- readH5AD(temp, use_hdf5 = TRUE)
    expect_s4_class(DelayedArray::seed(assay(backed)), "HDF5ArraySeed")
    expect_s4_class(assay(backed, "more"), "CsparseMatrix")
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

test_that("Reading H5AD works with native reader", {
    sce <- readH5AD(file, reader = "R")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading v0.8 H5AD works with native reader", {
    sce_py <- readH5AD(file_v08)
    sce_r <- readH5AD(file_v08, reader="R")

    expect_identical(rownames(sce_py), rownames(sce_r))
    expect_identical(colnames(sce_py), colnames(sce_r))

    expect_identical(rowData(sce_py), rowData(sce_r))

    expect_identical(colnames(colData(sce_py)), colnames(colData(sce_r)))
    expect_equal(colData(sce_py), colData(sce_r))

    # check the X assay
    expect_identical(assays(sce_py), assays(sce_r))

    # check the easy metadata columns
    for (key in c("dummy_category", "dummy_int", "dummy_int2", "highlight",
                  "iroot")) {
        expect_equal(metadata(sce_py)[[key]], metadata(sce_r)[[key]])
    }

    # For these columns the Python reader reads an array
    for (key in c("dummy_bool", "dummy_bool2")) {
        expect_equal(as.vector(metadata(sce_py)[[key]]), metadata(sce_r)[[key]])
    }
})

test_that("Skipping slot conversion works", {
    sce <- readH5AD(file, layers = FALSE, uns = FALSE, var = FALSE, obs = FALSE,
                    varm = FALSE, obsm = FALSE, varp = FALSE, obsp = FALSE)

    expect_identical(assayNames(sce), "X")
    expect_identical(metadata(sce), list())
    expect_equal(ncol(rowData(sce)), 0)
    expect_equal(ncol(colData(sce)), 0)
    expect_equal(length(reducedDims(sce)), 0)
    expect_equal(length(rowPairs(sce)), 0)
    expect_equal(length(colPairs(sce)), 0)
})

test_that("Selective slot conversion works", {
    sce <- readH5AD(file, uns = "iroot")

    expect_identical(names(metadata(sce)), "iroot")
})

test_that("Selective DF conversion works", {
    sce <- readH5AD(file, obs = "cell_type")

    expect_identical(names(colData(sce)), "cell_type")
})

test_that("Conversion of raw works", {
    skip_if_offline()

    cache <- BiocFileCache::BiocFileCache(ask = FALSE)
    example_file <- BiocFileCache::bfcrpath(
        cache, "https://ndownloader.figshare.com/files/30462915"
    )

    sce <- readH5AD(example_file, raw = TRUE)

    names <- list(
        assays = c("X"),
        colData = c("n_genes", "n_genes_by_counts", "total_counts",
                    "total_counts_mt", "pct_counts_mt", "leiden"),
        rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts",
                    "mean_counts", "pct_dropout_by_counts", "total_counts",
                    "highly_variable", "means", "dispersions",
                    "dispersions_norm", "mean", "std"),
        metadata = c("hvg", "leiden", "neighbors", "pca", "rank_genes_groups",
                     "umap"),
        redDim = c("X_pca", "X_umap"),
        varm = c("PCs"),
        colPairs = c("connectivities", "distances"),
        raw_rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts",
                        "mean_counts", "pct_dropout_by_counts", "total_counts",
                        "highly_variable", "means", "dispersions",
                        "dispersions_norm")
    )

    missing <- list()

    validateH5ADSCE(sce, names, missing)
})
