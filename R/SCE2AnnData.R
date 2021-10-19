#' @rdname AnnData-Conversion
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param X_name For `SCE2AnnData()` name of the assay to use as the primary
#' matrix (`X`) of the AnnData object. If `NULL`, the first assay of `sce` will
#' be used by default. For `AnnData2SCE()` name used when saving `X` as an
#' assay. If `NULL` looks for an `X_name` value in `uns`, otherwise uses `"X"`.
#' @param verbose Logical scalar indicating whether to print progress messages.
#' If `NULL` uses `getOption("zellkonverter.verbose")`.
#'
#' @export
#' @importFrom utils capture.output
#' @importFrom S4Vectors metadata
#' @importFrom reticulate import r_to_py
SCE2AnnData <- function(sce, X_name = NULL, skip_assays = FALSE,
                        verbose = NULL) {
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
    adata <- anndata$AnnData(X = X)
    cli::cli_progress_done()

    assay_names <- assayNames(sce)
    assay_names <- assay_names[!assay_names == X_name]
    if (length(assay_names) > 0) {
        .ui_step(
            "Converting {.field additional assays} to {.field layers}",
            msg_done = "{.field additional assays} converted to {.field layers}"
        )
        if (!skip_assays) {
            assays_list <- assays(sce, withDimnames = FALSE)
            assays_list <- lapply(assays_list[assay_names], .makeNumpyFriendly)
        } else {
            assays_list <- rep(list(fake_mat), length(assay_names))
            names(assays_list) <- assay_names
        }
        adata$layers <- assays_list
        cli::cli_progress_done()
    } else {
        .ui_info("No {.field additional assays} present")
    }

    col_data <- colData(sce)
    is_atomic <- vapply(col_data, is.atomic, NA)
    if (any(!is_atomic)) {
        non_atomic_cols <- colnames(col_data)[!is_atomic]
        warning(
            "The following colData columns are not atomic and will be stored ",
            "in metadata(sce)$.colData before conversion: ",
            paste(non_atomic_cols, collapse = ", ")
        )

        if (".colData" %in% names(metadata(sce))) {
            meta_list <- metadata(sce)$.colData
        } else {
            meta_list <- list()
        }

        for (col in non_atomic_cols) {
            store_name <- make.names(c(col, names(meta_list)), unique = TRUE)[1]
            meta_list[[store_name]] <- col_data[[col]]
        }

        col_data[non_atomic_cols] <- NULL
        metadata(sce)$.colData <- meta_list
    }

    if (ncol(col_data) > 0) {
        .ui_step(
            "Converting {.field colData} to {.field obs}",
            msg_done = "{.field colData} converted to {.field obs}"
        )
        # Manually construct the data.frame to avoid mangling column names
        obs <- do.call(
            data.frame,
            c(
                as.list(col_data),
                check.names      = FALSE,
                stringsAsFactors = FALSE
            )
        )
        adata$obs <- obs
        cli::cli_progress_done()
    } else {
        .ui_info("{.field colData} is empty and was skipped")
    }

    row_data <- rowData(sce)
    if (!is.null(int_metadata(sce)$has_varm)) {
        .ui_step(
            "Converting {.field rowData$varm} to {.field varm}",
            msg_done = "{.field rowData$varm} converted to {.field varm}"
        )
        varm <- as.list(row_data[["varm"]])
        row_data[["varm"]] <- NULL
        adata$varm <- varm
        cli::cli_progress_done()
    } else {
        .ui_info("{.field rowData$varm} is empty and was skipped")
    }

    is_atomic <- vapply(row_data, is.atomic, NA)
    if (any(!is_atomic)) {
        non_atomic_cols <- colnames(row_data)[!is_atomic]
        warning(
            "The following rowData columns are not atomic and will be stored ",
            "in metadata(sce)$.rowData before conversion: ",
            paste(non_atomic_cols, collapse = ", ")
        )

        if (".rowData" %in% names(metadata(sce))) {
            meta_list <- metadata(sce)$.rowData
        } else {
            meta_list <- list()
        }

        for (col in non_atomic_cols) {
            store_name <- make.names(c(col, names(meta_list)), unique = TRUE)[1]
            meta_list[[store_name]] <- row_data[[col]]
        }

        row_data[non_atomic_cols] <- NULL
        metadata(sce)$.rowData <- meta_list
    }

    if (ncol(row_data) > 0) {
        .ui_step(
            "Converting {.field rowData} to {.field var}",
            msg_done = "{.field rowData} converted to {.field var}"
        )
        # Manually construct the data.frame to avoid mangling column names
        var <- do.call(
            data.frame,
            c(
                as.list(row_data),
                check.names      = FALSE,
                stringsAsFactors = FALSE
            )
        )
        adata$var <- var
        cli::cli_progress_done()
    } else {
        .ui_info("{.field rowData} is empty and was skipped")
    }

    red_dims <- as.list(reducedDims(sce))
    if (length(red_dims) > 0) {
        .ui_step(
            "Converting {.field reducedDims} to {.field obsm}",
            msg_done = "{.field reducedDims} converted to {.field obsm}"
        )
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
    } else {
        .ui_info("{.field reducedDims} is empty and was skipped")
    }

    meta_list <- metadata(sce)
    if (length(meta_list) > 0) {
        .ui_step(
            "Converting {.field metadata} to {.field uns}",
            msg_done = "{.field metadata} converted to {.field uns}"
        )
        meta_list <- .addListNames(meta_list)
        uns_list <- list()
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
        uns_list[["X_name"]] <- X_name
        adata$uns <- reticulate::dict(uns_list)
        cli::cli_progress_done()
    } else {
        .ui_info("{.field metadata} is empty and was skipped")
    }

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

    if (length(colPairs(sce)) > 0) {
        .ui_step(
            "Converting {.field colPairs} to {.field obsp}",
            msg_done = "{.field colPairs} converted to {.field obsp}"
        )
        adata$obsp <- as.list(colPairs(sce, asSparse = TRUE))
        cli::cli_progress_done()
    } else {
        .ui_info("{.field colPairs} is empty and was skipped")
    }

    if (!is.null(colnames(sce))) {
        adata$obs_names <- colnames(sce)
    }

    if (!is.null(rownames(sce))) {
        adata$var_names <- rownames(sce)
    }

    adata
}

#' @importFrom methods as is
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom DelayedArray is_sparse
#' @importFrom Matrix t
.makeNumpyFriendly <- function(x, transpose = TRUE) {
    if (transpose) {
        x <- t(x)
    }

    # Code from Charlotte Soneson in kevinrue/velociraptor.
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
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
