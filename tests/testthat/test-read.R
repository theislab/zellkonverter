test_that("Reading H5AD works", {
    file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")
})
