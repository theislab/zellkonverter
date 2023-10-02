# This file tests compatibility with the R {anndata} package
# Despite best efforts the package isn't reliably unloaded so these tests have
# been moved to a separate file that is (hopefully) always run last

test_that("Reading is compatible with R anndata", {
    skip_if_offline()
    skip_if_not_installed("withr")
    skip_if_not_installed("anndata")

    withr::with_package("anndata", {
        file <- system.file("extdata", "krumsiek11.h5ad",
                            package = "zellkonverter")
        sce <- readH5AD(file)
        expect_s4_class(sce, "SingleCellExperiment")

        expect_identical(assayNames(sce), "X")
        expect_identical(colnames(colData(sce)), "cell_type")

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
            raw_rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts",
                            "mean_counts", "pct_dropout_by_counts",
                            "total_counts", "highly_variable", "means",
                            "dispersions", "dispersions_norm"),
            redDim = c("X_pca", "X_umap"),
            varm = c("PCs"),
            colPairs = c("connectivities", "distances"),
            metadata = c("hvg", "leiden", "neighbors", "pca",
                         "rank_genes_groups", "umap")
        )

        missing <- list()

        validateH5ADSCE(sce, names, missing)
    })

    pkgload::unload("anndata")
})

test_that("Writing is compatible with R anndata", {
    skip_if_offline()
    skip_if_not_installed("withr")
    skip_if_not_installed("anndata")

    withr::with_package("anndata", {
        sce <- scRNAseq::ZeiselBrainData()
        temp <- tempfile(fileext = ".h5ad")
        writeH5AD(sce, temp)
        expect_true(file.exists(temp))

        # Reading it back out again. Hopefully we didn't lose anything important
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

    pkgload::unload("anndata")
})
