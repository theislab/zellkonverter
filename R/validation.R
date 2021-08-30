#' Validate H5AD SCE
#'
#' Validate a SingleCellExperiment created by `readH5AD()`. Designed to be used
#' inside `testhat::test_that()` during package testing.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param names Named list of expected names. Names are slots and values are
#' vectors of names that are expected to exist in that slot.
#' @param missing Named list of known missing names. Names are slots and values
#' are vectors of names that are expected to not exist in that slot.
#'
#' @details
#' This function checks that a SingleCellExperiment contains the expected items
#' in each slot. The main reason for this function is avoid repeating code when
#' testing multiple `.h5ad` files. The following items in `names` and `missing`
#' are recognised:
#'
#' * `assays` - Assay names
#' * `colData` - colData column names
#' * `rowData` - rowData column names
#' * `metadata` - metadata names
#' * `redDim` - Reduced dimension names
#' * `varm` - Column names of the `varm` rowData column (from the AnnData varm
#'    slot)
#' * `colPairs` - Column pair names
#' * `rowPairs` - rowData pair names
#'
#' If an item in `names` or `missing` is `NULL` then it won't be checked. The
#' items in `missing` are checked that they explicitly do not exist. This is
#' mostly for record keeping when something is known to not be converted but can
#' also be useful when the corresponding `names` item is `NULL`.
#'
#' @return If checks are successful `TRUE` invisibly, if not other output
#' depending on the context
#'
#' @author Luke Zappia
#'
#' @examples
#' file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
#' sce <- readH5AD(file)
#'
#' names <- list(
#'     assays = "X",
#'     colData = "cell_type",
#'     metadata = c("highlights", "iroot"),
#' )
#' missing <- list(metadata = "something missing")
#'
#' validateH5ADSCE(sce, names, missing)
validateH5ADSCE <- function(sce, names, missing) {

    if ("varm" %in% colnames(SummarizedExperiment::rowData(sce))) {
        varm <- SummarizedExperiment::rowData(sce)$varm
        SummarizedExperiment::rowData(sce)$varm <- NULL
    } else {
        varm <- NULL
    }

    .names_validator(
        "Assay names",
        SummarizedExperiment::assayNames(sce),
        names$assays,
        missing$assays
    )

    .names_validator(
        "colData names",
        colnames(SummarizedExperiment::colData(sce)),
        names$colData,
        missing$colData
    )

    .names_validator(
        "rowData names",
        colnames(SummarizedExperiment::rowData(sce)),
        names$rowData,
        missing$rowData
    )

    .names_validator(
        "metadata names",
        names(S4Vectors::metadata(sce)),
        names$metadata,
        missing$metadata
    )

    .names_validator(
        "redDim names",
        SingleCellExperiment::reducedDimNames(sce),
        names$redDim,
        missing$redDim
    )

    .names_validator(
        "varm names",
        colnames(varm),
        names$varm,
        missing$varm
    )

    .names_validator(
        "colPairs names",
        names(SingleCellExperiment::colPairs(sce)),
        names$colPairs,
        missing$colPairs
    )

    .names_validator(
        "rowPairs names",
        names(SingleCellExperiment::rowPairs(sce)),
        names$rowPairs,
        missing$rowPairs
    )

    invisible(TRUE)
}

.names_validator <- function(label, actual_names, correct_names,
                             missing_names) {

    if (!is.null(correct_names)) {
        testthat::expect_identical(
            actual_names,
            correct_names,
            label = label
        )
    }

    if (!is.null(missing_names)) {
        testthat::expect_true(
            !any(missing_names %in% actual_names),
            label = paste(label, "missing")
        )
    }

    invisible(TRUE)
}

#' Expect SCE
#'
#' Test that a SingleCellExperiment matches an expected object. Designed to be
#' used inside `testhat::test_that()` during package testing.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param expected A template \linkS4class{SingleCellExperiment} object to
#' compare to.
#'
#' @return `TRUE` invisibly if checks pass
#'
#' @author Luke Zappia
#'
#' @examples
#' file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
#' sce1 <- readH5AD(file)
#' sce2 <- readH5AD(file)
#'
#' expectSCE(sce1, sce2))
expectSCE <- function(sce, expected) {
    expect_identical(dimnames(sce), dimnames(expected))
    if (length(metadata(expected)) > 0) {
        expect_identical(metadata(sce), metadata(expected))
    }
    expect_identical(assayNames(sce), assayNames(expected))
    for (assay in assayNames(expected)) {
        expect_equal(assay(sce, assay), assay(expected, assay))
    }
    expect_identical(reducedDims(sce), reducedDims(expected))
    expect_identical(rowData(sce), rowData(expected))
    expect_identical(colData(sce), colData(expected))
    expect_identical(rowPairs(sce), rowPairs(expected))
    expect_identical(colPairs(sce), colPairs(expected))

    invisible(TRUE)
}
