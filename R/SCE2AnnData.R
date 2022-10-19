#' @rdname AnnData-Conversion
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name For `SCE2AnnData()` name of the assay to use as the primary
#' matrix (`X`) of the AnnData object. If `NULL`, the first assay of `sce` will
#' be used by default. For `AnnData2SCE()` name used when saving `X` as an
#' assay. If `NULL` looks for an `X_name` value in `uns`, otherwise uses `"X"`.
#' @param assays,colData,rowData,reducedDims,metadata,colPairs,rowPairs
#' Arguments specifying how these slots are converted. If `TRUE` everything in
#' that slot is converted, if `FALSE` nothing is converted and if a character
#' vector only those items or columns are converted.
#' @param verbose Logical scalar indicating whether to print progress messages.
#' If `NULL` uses `getOption("zellkonverter.verbose")`.
#'
#' @export
#' @importFrom utils capture.output
#' @importFrom S4Vectors metadata
#' @importFrom reticulate import r_to_py py_to_r
SCE2AnnData <- function(sce, X_name = NULL, assays = TRUE, colData = TRUE,
                        rowData = TRUE, varm = TRUE, reducedDims = TRUE,
                        metadata = TRUE, colPairs = TRUE, rowPairs = TRUE,
                        skip_assays = FALSE, verbose = NULL) {
    anndata <- import("anndata")

    .ui_process(
        "Converting {.field AnnData} to {.field SingleCellExperiment}"
    )

    if (is.null(X_name)) {
        .ui_step(
            "Selecting {.field X matrix}",
            msg_done = "Selected {.field X matrix}"
        )
        if (length(assays(sce)) == 0) {
            stop("'sce' does not contain any assays")
        }
        X_name <- assayNames(sce)[1]
        cli::cli_alert_info(
            "Using the {.field '{X_name}'} assay as the {.field X matrix}"
        )
        cli::cli_progress_done()
    }

    .ui_step(
        "Converting {.field assays${X_name}} to {.field X matrix}",
        msg_done = "{.field assays${X_name}} converted to {.field X matrix}"
    )
    if (!skip_assays) {
        X <- assay(sce, X_name)
        X <- .makeNumpyFriendly(X)
    } else {
        cli::cli_alert_warning(paste(
            "{.field skip_assays} is {.field TRUE}",
            "so {.field X/layers} will be empty"
        ))
        X <- fake_mat <- .make_fake_mat(rev(dim(sce)))
    }
    X <- reticulate::r_to_py(X)
    adata <- list(X = X, dtype = X$dtype)
    cli::cli_progress_done()

    assay_names <- assayNames(sce)
    assay_names <- assay_names[!assay_names == X_name]
    if (isFALSE(assays)) {
        .ui_info("Skipping conversion of {.field assays}")
    } else if (length(assay_names) == 0) {
        .ui_info("No {.field additional assays} present, assays were skipped")
    } else {
        .ui_step(
            "Converting {.field additional assays} to {.field layers}",
            msg_done = "{.field additional assays} converted to {.field layers}"
        )
        if (is.character(assays)) {
            if (!all(assays %in% assay_names)) {
                missing <- assays[!c(assays %in% assay_names)]
                .ui_warn(
                    "These selected assays are not in the object: {.field {missing}}"
                )
                warning(
                    "These selected assays are not in the object: ",
                    paste(missing, collapse = ", "),
                    call. = FALSE
                )
            }
            assay_names <- assay_names[assay_names %in% assays]
        }
        if (!skip_assays) {
            assays_list <- assays(sce, withDimnames = FALSE)
            assays_list <- lapply(assays_list[assay_names], .makeNumpyFriendly)
        } else {
            assays_list <- rep(list(fake_mat), length(assay_names))
            names(assays_list) <- assay_names
        }
        adata$layers <- assays_list
        cli::cli_progress_done()
    }

    if (isFALSE(colData)) {
        .ui_info("Skipping conversion of {.field colData}")
    } else {
        sce <- .store_non_atomic(sce, "colData")
        obs <- .convert_sce_df(colData(sce), "colData", "obs", select = colData)
        if (!is.null(obs)) {
            adata$obs <- obs
        }
    }

    if (!is.null(int_metadata(sce)$has_varm)) {
        varm_list <- as.list(rowData(sce)[["varm"]])
        rowData(sce)[["varm"]] <- NULL

        if (isFALSE(varm)) {
            .ui_info("Skipping conversion of {.field rowData$varm}")
        } else {
            .ui_step(
                "Converting {.field rowData$varm} to {.field varm}",
                msg_done = "{.field rowData$varm} converted to {.field varm}"
            )

            if (is.character(varm)) {
                varm <- .check_select(varm, "rowData$varm", names(varm_list))
                varm_list <- varm_list[varm]
            }

            # Make sure var names are set in case varm contains a data.frame
            if (!is.null(rownames(sce))) {
                adata$var_names <- rownames(sce)
            }

            adata$varm <- varm_list
            cli::cli_progress_done()
        }

    } else {
        .ui_info("{.field rowData$varm} is empty and was skipped")
    }

    if (isFALSE(rowData)) {
        .ui_info("Skipping conversion of {.field rowData}")
    } else {
        sce <- .store_non_atomic(sce, "rowData")
        var <- .convert_sce_df(rowData(sce), "rowData", "var", select = rowData)
        if (!is.null(var)) {
            adata$var <- var
        }
    }

    if (isFALSE(reducedDims)) {
        .ui_info("Skipping conversion of {.field reducedDims}")
    } else if (length(reducedDims(sce)) == 0) {
        .ui_info("{.field reducedDims} is empty and was skipped")
    } else {
        .ui_step(
            "Converting {.field reducedDims} to {.field obsm}",
            msg_done = "{.field reducedDims} converted to {.field obsm}"
        )
        red_dims <- as.list(reducedDims(sce))
        if (is.character(reducedDims)) {
            reducedDims <- .check_select(reducedDims, "reducedDims",
                                         names(red_dims))
            red_dims <- red_dims[reducedDims]
        }
        red_dims <- lapply(red_dims, .makeNumpyFriendly, transpose = FALSE)
        red_dims <- lapply(red_dims, function(rd) {
            if (!is.null(colnames(rd))) {
                rd <- r_to_py(as.data.frame(rd))
                rd <- rd$set_index(adata$obs_names)
            }

            rd
        })
        adata$obsm <- red_dims
        cli::cli_progress_done()
    }

    uns_list <- list()
    uns_list[["X_name"]] <- X_name
    if (isFALSE(metadata)) {
        .ui_info("Skipping conversion of {.field metadata}")
    } else if (length(metadata(sce)) == 0) {
        .ui_info("{.field metadata} is empty and was skipped")
    } else {
        .ui_step(
            "Converting {.field metadata} to {.field uns}",
            msg_done = "{.field metadata} converted to {.field uns}"
        )
        meta_list <- .addListNames(metadata(sce))
        if (is.character(metadata)) {
            metadata <- .check_select(metadata, "metadata", names(meta_list))
            meta_list <- meta_list[metadata]
        }
        for (item_name in names(meta_list)) {
            item <- meta_list[[item_name]]
            tryCatch(
                {
                    # Try to convert the item using reticulate, skip if it fails
                    # Capture the object output printed by reticulate
                    capture.output(r_to_py(item))
                    uns_list[[item_name]] <- item
                },
                error = function(err) {
                    warning(
                        "the '", item_name, "' item in 'metadata' cannot be ",
                        "converted to a Python type and has been skipped"
                    )
                }
            )
        }
        cli::cli_progress_done()
    }
    adata$uns <- reticulate::dict(uns_list)

    if (length(rowPairs(sce)) > 0) {
        .ui_step(
            "Converting {.field rowPairs} to {.field varp}",
            msg_done = "{.field rowPairs} converted to {.field varp}"
        )
        adata$varp <- as.list(rowPairs(sce, asSparse = TRUE))
        cli::cli_progress_done()
    } else {
        .ui_info("{.field rowPairs} is empty and was skipped")
    }

    adata$obsp <- .convert_sce_pairs(sce, "colPairs", "obsp", colPairs)

    adata$varp <- .convert_sce_pairs(sce, "rowPairs", "varp", rowPairs)

    if (!is.null(colnames(sce))) {
        if (is.null(adata$obs)) {
            adata$obs <- data.frame(row.names = colnames(sce))
        } else {
            rownames(adata$obs) <- colnames(sce)
        }
    }

    if (!is.null(rownames(sce))) {
        if (is.null(adata$var)) {
            adata$var <- data.frame(row.names = rownames(sce))
        } else {
            rownames(adata$var) <- rownames(sce)
        }
    }

    do.call(anndata$AnnData, adata)
}

#' @importFrom methods as is
#' @importClassesFrom Matrix CsparseMatrix
#' @importFrom DelayedArray is_sparse
#' @importFrom Matrix t
.makeNumpyFriendly <- function(x, transpose = TRUE) {
    if (transpose) {
        x <- t(x)
    }

    # Code from Charlotte Soneson in kevinrue/velociraptor.
    if (is_sparse(x)) {
        as(x, "CsparseMatrix")
    } else {
        as.matrix(x)
    }
}

.addListNames <- function(x) {
    if (length(x) == 0) {
        return(x)
    }

    if (is.null(names(x))) {
        names(x) <- paste0("item", seq_along(x))
        return(x)
    }

    list_names <- names(x)
    is_empty <- list_names == ""
    list_names[is_empty] <- paste0("item", seq_along(x)[is_empty])
    list_names <- make.names(list_names, unique = TRUE)

    names(x) <- list_names

    return(x)
}

.store_non_atomic <- function(sce, slot = c("rowData", "colData")) {

    slot <- match.arg(slot)

    df <- switch (slot,
        rowData = rowData(sce),
        colData = colData(sce)
    )

    is_atomic <- vapply(df, is.atomic, NA)

    if (all(is_atomic)) {
        return(sce)
    }

    non_atomic_cols <- colnames(df)[!is_atomic]
    warning(
        "The following ", slot, " columns are not atomic and will be ",
        "stored in metadata(sce)$.colData before conversion: ",
        paste(non_atomic_cols, collapse = ", ")
    )

    meta_slot <- paste0(".", slot)
    if (meta_slot %in% names(metadata(sce))) {
        meta_list <- metadata(sce)[[meta_slot]]
    } else {
        meta_list <- list()
    }

    for (col in non_atomic_cols) {
        store_name <- make.names(c(col, names(meta_list)), unique = TRUE)[1]
        meta_list[[store_name]] <- df[[col]]
    }

    df[non_atomic_cols] <- NULL
    metadata(sce)[[meta_slot]] <- meta_list

    if (slot == "rowData") {
        rowData(sce) <- df
    } else (
        colData(sce) <- df
    )

    return(sce)
}

.check_select <- function(select, slot_name, options) {

    verbose <- parent.frame()$verbose

    if (!all(select %in% options)) {
        missing <- select[!c(select %in% options)]
        .ui_warn(paste(
            "These selected {.field {slot_name}} items are not in the",
            "object: {.field {missing}}"
        ))
        warning(
            "These selected ", slot_name, " items are not in the ",
            "object: ", paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    select <- select[select %in% options]

    return(select)
}

.convert_sce_df <- function(sce_df, slot_name, to_name, select = TRUE) {

    if (ncol(sce_df) == 0) {
        .ui_info("{.field {slot_name}} is empty and was skipped")
        return(NULL)
    }

    .ui_step(
        "Converting {.field {slot_name}} to {.field {to_name}}",
        msg_done = "{.field {slot_name}} converted to {.field {to_name}}"
    )
    if (is.character(select)) {

        select <- .check_select(select, slot_name, colnames(sce_df))

        if (length(select) == 0) {
            return(NULL)
        }

        df <- sce_df[, select, drop = FALSE]
    } else {
        df <- sce_df
    }

    df <- do.call(
        data.frame,
        c(
            as.list(df),
            check.names      = FALSE,
            stringsAsFactors = FALSE
        )
    )
    cli::cli_progress_done()

    return(df)
}

.convert_sce_pairs <- function(sce, slot_name = c("rowPairs", "colPairs"),
                               to_name, select) {

    slot_name <- match.arg(slot_name)


    if (isFALSE(select)) {
        .ui_info("Skipping conversion of {.field {slot_name}}")
        return(list())
    }

    pairs <- switch (slot_name,
        rowPairs = as.list(rowPairs(sce, asSparse = TRUE)),
        colPairs = as.list(colPairs(sce, asSparse = TRUE))
    )

    if (length(pairs) == 0) {
        .ui_info("{.field {slot_name}} is empty and was skipped")
        return(list())
    }

    .ui_step(
        "Converting {.field {slot_name}} to {.field {to_name}}",
        msg_done = "{.field {slot_name}} converted to {.field {to_name}}"
    )

    if (is.character(select)) {
        select <- .check_select(select, slot_name, names(pairs))
        pairs <- pairs[select]
    }
    cli::cli_progress_done()

    return(pairs)
}
