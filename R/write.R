#' Write H5AD
#'
#' Write a H5AD file from a \linkS4class{SingleCellExperiment} object.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param file String containing a path to write the new `.h5ad` file.
#' @param X_name Name of the assay to use as the primary matrix (`X`) of the
#' AnnData object. If `NULL`, the first assay of `sce` will be used by default.
#' @param skip_assays Logical scalar indicating whether assay matrices should
#' be ignored when writing to `file`.
#' @param compression Type of compression when writing the new `.h5ad` file.
#' @param version A string giving the version of the **anndata** Python library
#' to use. Allowed values are available in `.AnnDataVersions`. By default the
#' latest version is used.
#' @param verbose Logical scalar indicating whether to print progress messages.
#' If `NULL` uses `getOption("zellkonverter.verbose")`.
#' @inheritDotParams SCE2AnnData
#'
#' @details
#'
#' ## Skipping assays
#'
#' Setting `skip_assays = TRUE` can occasionally be useful if the matrices in
#' `sce` are stored in a format that is not amenable for efficient conversion
#' to a **numpy**-compatible format. In such cases, it can be better to create
#' an empty placeholder dataset in `file` and fill it in R afterwards.
#'
#' ## **DelayedArray** assays
#'
#' If `sce` contains any **DelayedArray** matrices as assays `writeH5AD()` will
#' write them to disk using the **rhdf5** package directly rather than via
#' Python to avoid instantiating them in memory. However there is currently
#' an issue which prevents this being done for sparse **DelayedArray** matrices.
#'
#' ## Known conversion issues
#'
#' ### Coercion to factors
#'
#' The **anndata** package automatically converts some character vectors to
#' factors when saving `.h5ad` files. This can effect columns of `rowData(sce)`
#' and `colData(sce)` which may change type when the `.h5ad` file is read back
#' into R.
#'
#' ## Environment
#'
#' See [AnnData-Environment] for more details on **zellkonverter** Python
#' environments.
#'
#' @return A `NULL` is invisibly returned.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @seealso
#' [`readH5AD()`], to read a \linkS4class{SingleCellExperiment} file from a H5AD
#' file.
#'
#' [`SCE2AnnData()`], for developers to create an AnnData object from a
#' \linkS4class{SingleCellExperiment}.
#'
#' @examples
#' # Using the Zeisel brain dataset
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(scRNAseq)
#'     sce <- ZeiselBrainData()
#'
#'     # Writing to a H5AD file
#'     temp <- tempfile(fileext = ".h5ad")
#'     writeH5AD(sce, temp)
#' }
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom Matrix sparseMatrix
#' @importFrom DelayedArray is_sparse
writeH5AD <- function(sce, file, X_name = NULL, skip_assays = FALSE,
                      compression = c("none", "gzip", "lzf"), version = NULL,
                      verbose = NULL, ...) {

    compression <- match.arg(compression)

    if (compression == "none") {
        compression <- NULL
    }

    # Loop over and replace DelayedArrays.
    ass_list <- assays(sce)
    is_da <- logical(length(ass_list))
    for (a in seq_along(ass_list)) {
        # Skip sparse DelayedArrays due to rhdf5 issue
        # https://github.com/grimbough/rhdf5/issues/79
        if (is(ass_list[[a]], "DelayedMatrix") && !is_sparse(ass_list[[a]])) {
            is_da[a] <- TRUE
            assay(sce, a, withDimnames = FALSE) <- .make_fake_mat(dim(sce))
        }
    }

    env <- zellkonverterAnnDataEnv(version)
    version <- gsub("zellkonverterAnnDataEnv-", "", slot(env, "envname"))
    .ui_info("Using {.field anndata} version {.field {version}}")

    file <- path.expand(file)
    basiliskRun(
        env = env,
        fun = .H5ADwriter,
        sce = sce,
        file = file,
        X_name = X_name,
        skip_assays = skip_assays,
        compression = compression,
        verbose = verbose,
        ...
    )

    # Going back out and replacing each of them.
    if (any(is_da)) {
        for (p in which(is_da)) {
            if (p == 1L) {
                curp <- "X"
            } else {
                curp <- file.path("layers", assayNames(sce)[p])
            }
            rhdf5::h5delete(file, curp)
            mat <- ass_list[[p]]

            if (!is_sparse(mat)) {
                HDF5Array::writeHDF5Array(mat, filepath = file, name = curp)
            } else {
                .write_CSR_matrix(file, name = curp, mat = mat)
            }
        }
    }

    invisible(NULL)
}

#' @importFrom reticulate import
.H5ADwriter <- function(sce, file, X_name, skip_assays, compression,
                        verbose = NULL, ...) {
    adata <- SCE2AnnData(sce, X_name = X_name, skip_assays = skip_assays,
                         verbose = verbose, ...)
    .ui_step(
        "Writing {.file { .trim_path(file)} }",
        msg_done = "Wrote {.file { .trim_path(file)} }",
        spinner = TRUE
    )
    if (!is.null(compression)) {
        .ui_info("Using {.field compression} compression")
    }
    adata$write_h5ad(file, compression = compression)
}

# nocov start

# Skipping code coverage on these function because they aren't used until the
# sparse DelayedArray rhdf5 issue mentioned above is addressed

#' @importFrom DelayedArray blockApply rowAutoGrid type
.write_CSR_matrix <- function(file, name, mat, chunk_dim = 10000) {
    handle <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(handle))

    rhdf5::h5createGroup(handle, name)
    ghandle <- rhdf5::H5Gopen(handle, name)
    on.exit(rhdf5::H5Gclose(ghandle), add = TRUE, after = FALSE)

    rhdf5::h5writeAttribute("csc_matrix", ghandle, "encoding-type")
    rhdf5::h5writeAttribute("0.1.0", ghandle, "encoding-version")
    rhdf5::h5writeAttribute(rev(dim(mat)), ghandle, "shape")

    rhdf5::h5createDataset(
        handle,
        file.path(name, "data"),
        dims = 0,
        maxdims = rhdf5::H5Sunlimited(),
        H5type = if (type(mat) == "integer") {
            "H5T_NATIVE_INT32"
        } else {
            "H5T_NATIVE_DOUBLE"
        },
        chunk = chunk_dim
    )

    rhdf5::h5createDataset(
        handle,
        file.path(name, "indices"),
        dims    = 0,
        maxdims = rhdf5::H5Sunlimited(),
        H5type  = "H5T_NATIVE_UINT32",
        chunk   = chunk_dim
    )

    env <- new.env() # persist the 'last' counter.
    env$last <- 0L
    out <- blockApply(
        mat,
        grid      = rowAutoGrid(mat),
        FUN       = .blockwise_sparse_writer,
        env       = env,
        file      = handle,
        name      = name,
        as.sparse = TRUE
    )

    out <- as.double(unlist(out))
    iname <- file.path(name, "indptr")

    rhdf5::h5createDataset(
        handle,
        iname,
        dims   = length(out) + 1L,
        H5type = "H5T_NATIVE_UINT64"
    )

    rhdf5::h5writeDataset(c(0, cumsum(out)), handle, iname)
}

#' @importFrom DelayedArray nzdata nzindex
.blockwise_sparse_writer <- function(block, env, file, name) {
    nzdex <- nzindex(block)
    i <- nzdex[, 1]
    j <- nzdex[, 2]
    v <- nzdata(block)

    o <- order(i)
    i <- i[o]
    j <- j[o]
    v <- v[o]

    last <- env$last
    index <- list(last + seq_along(j))

    iname <- file.path(name, "indices")
    rhdf5::h5set_extent(file, iname, last + length(j))
    rhdf5::h5writeDataset(j - 1L, file, iname, index = index)

    vname <- file.path(name, "data")
    rhdf5::h5set_extent(file, vname, last + length(j))
    rhdf5::h5writeDataset(v, file, vname, index = index)

    env$last <- last + length(j)
    tabulate(i, nrow(block))
}

# nocov end
