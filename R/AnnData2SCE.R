#' Convert between AnnData and SingleCellExperiment
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
#'     }, env = zellkonverterAnnDataEnv(), sce = seger)
#' }
#' @name AnnData-Conversion
#' @rdname AnnData-Conversion
NULL

#' @rdname AnnData-Conversion
#'
#' @param adata A **reticulate** reference to a Python AnnData object.
#' @param layers,uns,var,obs,varm,obsm,varp,obsp,raw Arguments specifying how
#' these slots are converted. If `TRUE` everything in that slot is converted, if
#' `FALSE` nothing is converted and if a character vector only those items or
#' columns are converted.
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
AnnData2SCE <- function(adata, X_name = NULL, layers = TRUE, uns = TRUE,
                        var = TRUE, obs = TRUE, varm = TRUE, obsm = TRUE,
                        varp = TRUE, obsp = TRUE, raw = FALSE,
                        skip_assays = FALSE, hdf5_backed = TRUE,
                        verbose = NULL) {

    # In case the user accidentally passes an AnnDataR6 object
    if (is(adata, "AnnDataR6")) {
        .ui_warn(paste(
            "The passed object is a 'AnnDataR6' object, conversion is likely ",
            "to be less reliable"
        ))
        adata <- r_to_py(adata)
    }

    # Disable automatic {reticulate} conversion for this object
    disable_conversion_scope(adata)

    .ui_process(
        "Converting {.field AnnData} to {.field SingleCellExperiment}"
    )

    py_builtins <- import_builtins()

    dims <- unlist(py_to_r(adata$shape))
    dims <- rev(dims)

    meta_list <- .convert_anndata_slot(
        adata,
        slot_name = "uns",
        slot_keys = py_builtins$list(adata$uns$keys()),
        to_name = "metadata",
        select = uns
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
        # do not apply py_to_r yet, because this is taken care of by
        # .extract_or_skip_assay(...)!
        mat = adata$X,
        name = "'X' matrix"
    )

    x_mat <- x_out$mat
    obs_names <- py_to_r(adata$obs_names$to_list())
    var_names <- py_to_r(adata$var_names$to_list())
    # DelayedArray won't accept an empty vector for dimnames so set to NULL
    if (length(obs_names) == 0) {
        obs_names <- NULL
    }
    if (length(var_names) == 0) {
        var_names <- NULL
    }
    colnames(x_mat) <- obs_names
    rownames(x_mat) <- var_names
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
    if (isFALSE(layers)) {
        .ui_info("Skipping conversion of {.field layers}")
    } else if (length(layer_names) == 0) {
        .ui_info("{.field layers} is empty and was skipped")
    } else {
        .ui_process("Converting {.field layers} to {.field assays}")
        if (is.character(layers)) {
            if (!all(layers %in% layer_names)) {
                missing <- layers[!c(layers %in% layer_names)]
                .ui_warn(
                    "These selected layers are not in the object: {.field {missing}}"
                )
            }
            layer_names <- layer_names[layer_names %in% layers]
        }
        for (layer_name in layer_names) {
            .ui_step(
                "Converting {.field layers${layer_name}}",
                msg_done = "{.field layers${layer_name}} converted"
            )
            layer_out <- .extract_or_skip_assay(
                skip_assays = skip_assays,
                hdf5_backed = hdf5_backed,
                dims = dims,
                # do not apply py_to_r yet, because this is taken care of by
                # .extract_or_skip_assay(...)!
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

    row_data <- .convert_anndata_df(
        py_to_r(adata$var),
        slot_name = "var",
        to_name   = "rowData",
        select    = var
    )

    col_data <- .convert_anndata_df(
        py_to_r(adata$obs),
        slot_name = "obs",
        to_name   = "colData",
        select    = obs
    )

    varm_list <- .convert_anndata_slot(
        adata,
        slot_name = "varm",
        slot_keys = py_to_r(adata$varm_keys()),
        to_name   = "rowData$varm",
        select    = varm
    )

    if (length(varm_list) > 0) {
        # Create an empty DataFrame with the correct number of rows
        varm_df <- make_zero_col_DFrame(py_to_r(adata$n_vars))
        for (varm_name in names(varm_list)) {
            varm_df[[varm_name]] <- varm_list[[varm_name]]
        }
        row_data$varm <- varm_df
    }

    reddim_list <- .convert_anndata_slot(
        adata, "obsm", py_to_r(adata$obsm_keys()), "reducedDims", select = obsm
    )
    reddim_list <- lapply(reddim_list, as.matrix)

    varp_list <- .convert_anndata_slot(
        adata, "varp", py_builtins$list(adata$varp$keys()), "rowPairs",
        select = varp
    )

    obsp_list <- .convert_anndata_slot(
        adata, "obsp", py_builtins$list(adata$obsp$keys()), "colPairs",
        select = obsp
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

    if (isFALSE(raw)) {
        .ui_info("Skipping conversion of {.field raw}")
    } else if (is.null(py_to_r(adata$raw))) {
        .ui_info("{.field raw} is empty and was skipped")
    } else {
        .ui_process("Converting {.field raw} to {.field altExp}")

        raw_x <- .extract_or_skip_assay(
            skip_assays = skip_assays,
            hdf5_backed = hdf5_backed,
            dims = rev(unlist(py_to_r(adata$raw$shape))),
            # do not apply py_to_r yet, because this is taken care of by
            # .extract_or_skip_assay(...)!
            mat = adata$raw$X,
            name = "raw 'X' matrix"
        )
        colnames(raw_x$mat) <- colnames(output)

        raw_rowData <- .convert_anndata_df(py_to_r(adata$raw$var), "raw var",
                                           "raw rowData", select = TRUE)

        raw_varm_list <- .convert_anndata_slot(
            adata, "varm", py_builtins$list(adata$raw$varm$keys()),
            "raw rowData$varm", select = TRUE, raw = TRUE
        )

        if (length(raw_varm_list) > 0) {
            # Create an empty DataFrame with the correct number of rows
            raw_varm_df <- make_zero_col_DFrame(nrow(raw_x))
            for (varm_name in names(raw_varm_list)) {
                raw_varm_df[[varm_name]] <- raw_varm_list[[varm_name]]
            }
            raw_rowData$varm <- raw_varm_df
        }

        altExp(output, "raw") <- SummarizedExperiment(
            assays = list(X = raw_x$mat),
            rowData = raw_rowData
        )

        .ui_process_done()
    }

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
        if (hdf5_backed &&
            (is(mat, "h5py._hl.dataset.Dataset") ||
             is(mat, "anndata._core.sparse_dataset.SparseDataset"))) {
            file <- as.character(py_to_r(mat$file$id$name))
            name <- as.character(py_to_r(mat$name))
            if (.h5isgroup(file, name)) {
                mat <- HDF5Array::H5SparseMatrix(file, name)
            } else {
                mat <- HDF5Array::HDF5Array(file, name)
            }
        } else {
            mat <- try(t(py_to_r(mat)), silent = TRUE)
            if (is(mat, "try-error")) {
                if (isFALSE(skip_assays)) {
                    .ui_warn(paste(
                        "{.field {name}} does not support transposition and",
                        "has been skipped"
                    ))
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

.convert_anndata_slot <- function(adata, slot_name, slot_keys, to_name,
                                  select = TRUE, raw = FALSE) {

    verbose <- parent.frame()$verbose

    if (isFALSE(select)) {
        .ui_info("Skipping conversion of {.field {slot_name}}")
        return(list())
    }

    if (length(slot_keys) == 0) {
        .ui_info("{.field {slot_name}} is empty and was skipped")
        return(list())
    }

    .ui_process("Converting {.field {slot_name}} to {.field {to_name}}")
    .ui_step(
        "Converting {.field {slot_name}}",
        msg_done = "{.field {slot_name}} converted"
    )
    if (is.character(select)) {
        if (!all(select %in% slot_keys)) {
            missing <- select[!c(select %in% slot_keys)]
            .ui_warn(paste(
                "These selected {.field {slot_name}} items are not in the ",
                "object: {.field {missing}}"
            ))
        }
        slot_keys <- slot_keys[slot_keys %in% select]
    }

    if (raw) {
        adata <- adata$raw
    }

    converted <- .convert_anndata_list(
        adata[[slot_name]],
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

    # Check if list items can be accessed (can fail when {anndata} is loaded)
    # Attempt whole list conversion is accessing individual items fails
    adata_list <- tryCatch(
        {
            # Check if an item can be accessed, if yes return the list
            adata_list[[keys[1]]]
            adata_list
        },
        error = function(err) {
            # If not issue a warning and try to convert the whole list
            .ui_warn(paste(
                "Unable to access items in {.field {parent}}, attempting to ",
                "convert the whole list.\n",
                "Access error message: {.val {err$message}}"
            ))
            adata_list <- tryCatch(
                {
                    # If conversion is successful return the whole list
                    adata_list <- py_to_r(adata_list)
                    adata_list
                },
                error = function(err) {
                    # If whole list conversion fails issue a warning and return
                    # NULL
                    .ui_warn(paste(
                        "Whole list conversion failed for {.field {parent}}, ",
                        "this slot will not be converted.\n",
                        "Conversion error message: {.val {err$message}}"
                    ))
                    NULL
                }
            )
            adata_list
        }
    )

    # If items cannot be accessed return empty list
    if (is.null(adata_list)) {
        return(list())
    }

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

                if (is(item, "python.builtin.dict")) {
                    item <- .convert_anndata_list(
                        item, paste(parent, key, sep = "$")
                    )
                }

                if (is(item, "python.builtin.object")) {
                    item <- py_to_r(item)

                    if (inherits(item, "list")) {
                        item <- .convert_anndata_list(
                            item, paste(parent, key, sep = "$")
                        )
                    }
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
                .ui_warn(paste(
                    "Conversion failed for the item {.field {key}}",
                    "in {.field {parent}} with the following error and has",
                    "been skipped\n",
                    "Conversion error message: {.val {err$message}}"
                ))

                "failed"
            }
        )
        cli::cli_progress_done(result = status)
    }

    orig_names <- names(converted_list)
    names(converted_list) <- make.names(names(converted_list))
    is_modified <- names(converted_list) != orig_names
    if (any(is_modified)) {
        modifications <- paste0(
            "'", orig_names[is_modified], "'",
            " -> ",
            "'", names(converted_list)[is_modified], "'"
        )
        .ui_warn(paste(
            "The names of these selected {.field {parent}} items have",
            "been modified to match R conventions: {.field {modifications}}"
        ))
    }

    return(converted_list)
}

.convert_anndata_df <- function(adata_df, slot_name, to_name, select = TRUE) {

    verbose <- parent.frame()$verbose

    if (isFALSE(select)) {
        .ui_info("Skipping conversion of {.field {slot_name}}")
        return(make_zero_col_DFrame(nrow(adata_df)))
    }

    .ui_step(
        "Converting {.field {slot_name}} to {.field {to_name}}",
        msg_done = "{.field {slot_name}} converted to {.field {to_name}}"
    )
    if (is.character(select)) {
        if (!all(select %in% colnames(adata_df))) {
            missing <- select[!c(select %in% colnames(adata_df))]
            .ui_warn(paste(
                "These selected {.field {slot_name}} columns are not in the",
                "object: {.field {missing}}"
            ))
            select <- setdiff(select, missing)
        }
    } else {
       select <- colnames(adata_df)
    }

    df <- adata_df[, select, drop = FALSE]

    # Return early if there are no columns as the next steps can break nrow
    if (ncol(df) == 0) {
        cli::cli_progress_done()
        return(DataFrame(df))
    }

    # Second conversion by column for types the {reticulate} misses
    # (mostly Pandas arrays)
    # Also convert 1D arrays to vectors
    df <- lapply(df, function(col) {
        if (is(col, "python.builtin.object")) {
            col <- reticulate::py_to_r(col)
        }
        if (is(col, "array") && is.na(ncol(col))) {
            col <- as.vector(col)
        }
        col
    })

    df <- DataFrame(df)

    is_modified <- colnames(df) != select
    if (any(is_modified)) {
        modifications <- paste0(
            "'", select[is_modified], "'",
            " -> ",
            "'", colnames(df)[is_modified], "'"
        )
        .ui_warn(paste(
            "The names of these selected {.field {slot_name}} columns have",
            "been modified to match R conventions: {.field {modifications}}"
        ))
    }

    cli::cli_progress_done()

    return(df)
}
