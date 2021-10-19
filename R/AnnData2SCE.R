#' Convert AnnData between and SingleCellExperiment
#'
#' Conversion between Python AnnData objects and
#' \linkS4class{SingleCellExperiment} objects.
#'
#' @details
#' These functions assume that an appropriate Python environment has already
#' been loaded. As such, they are largely intended for developer use, most
#' typically inside a **basilisk** context.
#'
#' The conversion is not entirely lossless. The current mapping is shown below
#' (also at <https://tinyurl.com/AnnData2SCE>):
#'
#' \if{html}{
#'     \figure{AnnData2SCE.png}{options: width=800, alt="SCE-AnnData map"}
#' }
#' \if{latex}{\figure{AnnData2SCE.png}{options: width=5in}}
#'
#' In `SCE2AnnData()`, matrices are converted to a **numpy**-friendly format.
#' Sparse matrices are converted to \linkS4class{dgCMatrix} objects while all
#' other matrices are converted into ordinary matrices. If `skip_assays = TRUE`,
#' empty sparse matrices are created instead and the user is expected to fill in
#' the assays on the Python side.
#'
#' For `AnnData2SCE()`, a warning is raised if there is no corresponding R
#' format for a matrix in the AnnData object, and an empty sparse matrix is
#' created instead as a placeholder. If `skip_assays = NA`, no warning is
#' emitted but variables are created in the [`int_metadata()`] of the output to
#' specify which assays were skipped.
#'
#' If `skip_assays = TRUE`, empty sparse matrices are created for all assays,
#' regardless of whether they might be convertible to an R format or not.
#' In both cases, the user is expected to fill in the assays on the R side,
#' see [`readH5AD()`] for an example.
#'
#' We attempt to convert between items in the \linkS4class{SingleCellExperiment}
#' [`metadata()`] slot and the `AnnData` `uns` slot. If an item cannot be
#' converted a warning will be raised.
#'
#' Values stored in the `varm` slot of an `AnnData` object are stored in a
#' column of [`rowData()`] in a \linkS4class{SingleCellExperiment}
#' as a \linkS4class{DataFrame} of matrices. If this column is present an
#' attempt is made to transfer this information when converting from
#' \linkS4class{SingleCellExperiment} to `AnnData`.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @return `AnnData2SCE()` will return a \linkS4class{SingleCellExperiment}
#' containing the equivalent data from `adata`.
#'
#' `SCE2AnnData()` will return a **reticulate** reference to an AnnData object
#' containing the content of `sce`.
#'
#' @seealso
#' [`writeH5AD()`] and [`readH5AD()`] for dealing directly with H5AD files.
#'
#' @examples
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(basilisk)
#'     library(scRNAseq)
#'     seger <- SegerstolpePancreasData()
#'
#'     # These functions are designed to be run inside
#'     # a specified Python environment
#'     roundtrip <- basiliskRun(fun = function(sce) {
#'         # Convert SCE to AnnData:
#'         adata <- zellkonverter::SCE2AnnData(sce)
#'
#'         # Maybe do some work in Python on 'adata':
#'         # BLAH BLAH BLAH
#'
#'         # Convert back to an SCE:
#'         zellkonverter::AnnData2SCE(adata)
#'     }, env = zellkonverterAnnDataEnv, sce = seger)
#' }
#' @name AnnData-Conversion
#' @rdname AnnData-Conversion
NULL

#' @rdname AnnData-Conversion
#'
#' @param adata A **reticulate** reference to a Python AnnData object.
#' @param skip_assays Logical scalar indicating whether to skip conversion of
#' any assays in `sce` or `adata`, replacing them with empty sparse matrices
#' instead.
#' @param hdf5_backed Logical scalar indicating whether HDF5-backed matrices
#' in `adata` should be represented as HDF5Array objects. This assumes that
#' `adata` is created with `backed="r"`.
#' @param verbose Logical scalar indicating whether to print progress messages.
#' If `NULL` uses `getOption("zellkonverter.verbose")`.
#'
#' @export
#' @importFrom methods selectMethod is
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom reticulate import_builtins
AnnData2SCE <- function(adata, X_name = NULL, skip_assays = FALSE,
                        hdf5_backed = TRUE, verbose = NULL) {

    .ui_process(
        "Converting {.field AnnData} to {.field SingleCellExperiment}"
    )

    py_builtins <- import_builtins()

    dims <- unlist(adata$shape)
    dims <- rev(dims)

    meta_list <- .convert_anndata_slot(
        adata, "uns", py_builtins$list(adata$uns$keys()), "metadata"
    )

    .ui_step(
        "Converting {.field X matrix} to {.field assay}",
        msg_done = "{.field X matrix} converted to {.field assay}"
    )
    if (skip_assays) {
        cli::cli_alert_warning(
            "{.field skip_assays} is {.field TRUE} so assays will be empty"
        )
    }
    x_out <- .extract_or_skip_assay(
        skip_assays = skip_assays,
        hdf5_backed = hdf5_backed,
        dims = dims,
        mat = adata$X,
        name = "'X' matrix"
    )

    x_mat <- x_out$mat
    colnames(x_mat) <- adata$obs_names$to_list()
    rownames(x_mat) <- adata$var_names$to_list()
    skipped_x <- x_out$skipped

    if (is.null(X_name)) {
        if ("X_name" %in% names(meta_list)) {
            X_name <- meta_list[["X_name"]]
            cli::cli_alert_info("Using stored X_name value {.field '{X_name}'}")
            meta_list[["X_name"]] <- NULL
        } else {
            X_name <- "X"
        }
    }
    cli::cli_progress_done()

    assays_list <- list()
    assays_list[[X_name]] <- x_mat

    layer_names <- names(py_builtins$dict(adata$layers))
    skipped_layers <- character(0)
    if (length(layer_names) == 0) {
        .ui_info("{.field layers} is empty and was skipped")
    } else {
        .ui_process("Converting {.field layers} to {.field assays}")
        for (layer_name in layer_names) {
            .ui_step(
                "Converting {.field layers${layer_name}}",
                msg_done = "{.field layers${layer_name}} converted"
            )
            layer_out <- .extract_or_skip_assay(
                skip_assays = skip_assays,
                hdf5_backed = hdf5_backed,
                dims = dims,
                mat = adata$layers$get(layer_name),
                name = sprintf("'%s' layer matrix", layer_name)
            )
            if (layer_out$skipped) {
                skipped_layers <- c(skipped_layers, layer_name)
            }
            assays_list[[layer_name]] <- layer_out$mat
            cli::cli_progress_done()
        }
        .ui_process_done()
    }

    .ui_step(
        "Converting {.field var} to {.field rowData}",
        msg_done = "{.field var} converted to {.field rowData}"
    )
    row_data <- DataFrame(adata$var)
    cli::cli_progress_done()

    .ui_step(
        "Converting {.field obs} to {.field colData}",
        msg_done = "{.field obs} converted to {.field colData}"
    )
    col_data <- DataFrame(adata$obs)
    cli::cli_progress_done()

    varm_list <- .convert_anndata_slot(
        adata, "varm", adata$varm_keys(), "rowData$varm"
    )

    if (length(varm_list) > 0) {
        # Create an empty DataFrame with the correct number of rows
        varm_df <- make_zero_col_DFrame(adata$n_vars)
        for (varm_name in names(varm_list)) {
            varm_df[[varm_name]] <- varm_list[[varm_name]]
        }
        row_data$varm <- varm_df
    }

    reddim_list <- .convert_anndata_slot(
        adata, "obsm", adata$obsm_keys(), "reducedDims"
    )
    reddim_list <- lapply(reddim_list, as.matrix)

    varp_list <- .convert_anndata_slot(
        adata, "varp", py_builtins$list(adata$varp$keys()), "rowPairs"
    )

    obsp_list <- .convert_anndata_slot(
        adata, "obsp", py_builtins$list(adata$obsp$keys()), "colPairs"
    )

    .ui_step(
        "Constructing {.field SingleCellExperiment}",
        msg_done = "{.field SingleCellExperiment} constructed"
    )
    output <- SingleCellExperiment(
        assays = assays_list,
        rowData = row_data,
        colData = col_data,
        reducedDims = reddim_list,
        metadata = meta_list,
        rowPairs = varp_list,
        colPairs = obsp_list
    )

    # Specifying which assays got skipped, if the skipping was variable.
    if (is.na(skip_assays)) {
        int_metadata(output)$skipped_x <- skipped_x
        int_metadata(output)$skipped_layers <- skipped_layers
    }

    if (length(varm_list) > 0) {
        int_metadata(output)$has_varm <- names(varm_list)
    }
    cli::cli_progress_done()

    output
}

#' @importFrom Matrix t
.extract_or_skip_assay <- function(skip_assays, hdf5_backed, dims, mat, name) {
    skipped <- FALSE

    if (isTRUE(skip_assays)) {
        # Value of 'mat' is never used so the promise never evaluates; thus,
        # skip_assays=TRUE avoids any actual transfer of content from Python.
        mat <- .make_fake_mat(dims)
    } else {
        if (hdf5_backed && is(mat, "python.builtin.object")) {
            file <- as.character(mat$file$id$name)
            name <- as.character(mat$name)
            if (.h5isgroup(file, name)) {
                mat <- HDF5Array::H5SparseMatrix(file, name)
            } else {
                mat <- HDF5Array::HDF5Array(file, name)
            }
        } else {
            mat <- try(t(mat), silent = TRUE)
            if (is(mat, "try-error")) {
                if (isFALSE(skip_assays)) {
                    warning(
                        name,
                        " does not support transposition and has been skipped"
                    )
                }
                mat <- .make_fake_mat(dims)
                skipped <- TRUE
            }
        }
    }

    if (is(mat, "dgRMatrix")) {
        mat <- as(mat, "CsparseMatrix")
    }

    list(mat = mat, skipped = skipped)
}

#' @importFrom Matrix sparseMatrix
.make_fake_mat <- function(dims) {
    sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = dims
    )
}

# Borrowed from HDF4Array
# https://github.com/Bioconductor/HDF5Array/blob/fb015cf1c789bbb905a6bf8af2c7b50a24a60795/R/h5utils.R#L66-L75
.h5isgroup <- function(filepath, name) {
    fid <- rhdf5::H5Fopen(filepath, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(fid))
    gid <- try(rhdf5::H5Gopen(fid, name), silent = TRUE)
    ans <- !inherits(gid, "try-error")
    if (ans) {
        rhdf5::H5Gclose(gid)
    }

    return(ans)
}

.convert_anndata_slot <- function(adata, slot_name, slot_keys, to_name = "MISSING") {

    verbose <- parent.frame()$verbose

    if (length(slot_keys) == 0) {
        .ui_info("{.field {slot_name}} is empty and was skipped")
        return(list())
    }

    .ui_process("Converting {.field {slot_name}} to {.field {to_name}}")
    .ui_step(
        "Converting {.field {slot_name}}",
        msg_done = "{.field {slot_name}} converted"
    )
    converted <- .convert_anndata_list(
        adata[slot_name],
        parent = slot_name,
        keys = slot_keys
    )
    cli::cli_progress_done()

    return(converted)
}

.convert_anndata_list <- function(adata_list, parent,
                                  keys = names(adata_list)) {
    py_builtins <- import_builtins()

    converted_list <- list()

    verbose <- parent.frame()$verbose

    for (key in keys) {
        .ui_step(
            "Converting {.field {parent}${key}}",
            msg_done = "{.field {parent}${key}} converted"
        )
        status <- tryCatch(
            {
                item <- adata_list[[key]]

                item_type <- py_builtins$str(py_builtins$type(item))
                if (grepl("OverloadedDict", item_type)) {
                    item <- py_builtins$dict(item)
                }

                if (is(item, "python.builtin.object")) {
                    item <- reticulate::py_to_r(item)
                }

                if (inherits(item, "list")) {
                    item <- .convert_anndata_list(
                        item, paste(parent, key, sep = "$")
                    )
                }

                if (is.data.frame(item)) {
                    # Remove pandas index stored by reticulate which (should)
                    # be redundant with rownames as H5AD doesn't support
                    # multiple indexes (yet)
                    attr(item, "pandas.index") <- NULL
                }

                converted_list[[key]] <- item

                "done"
            },
            error = function(err) {
                warning(
                    "conversion failed for the item '",
                    key, "' in '", parent, "' with ",
                    "the following error and has been skipped\n",
                    "Conversion error message: ", err,
                    call. = FALSE
                )

                "failed"
            }
        )
        cli::cli_progress_done(result = status)
    }

    return(converted_list)
}
