#' Set zellkonverter verbose
#'
#' Set the zellkonverter verbosity option
#'
#' @param verbose Logical value for the verbosity option.
#'
#' @details
#' Running `setZellkonverterVerbose(TRUE)` will turn on **zellkonverter**
#' progress messages by default without having to set `verbose = TRUE` in each
#' function call. This is done by setting the `"zellkonverter.verbose"` option.
#' Running `setZellkonverterVerbose(FALSE)` will turn default verbosity off.
#'
#' @return The value of getOption("zellkonverter.verbose") invisibly
#' @export
#'
#' @examples
#' current <- getOption("zellkonverter.verbose")
#' setZellkonverterVerbose(TRUE)
#' getOption("zellkonverter.verbose")
#' setZellkonverterVerbose(FALSE)
#' getOption("zellkonverter.verbose")
#' setZellkonverterVerbose(current)
#' getOption("zellkonverter.verbose")
setZellkonverterVerbose <- function(verbose = TRUE) {
    options(zellkonverter.verbose = isTRUE(verbose))
    invisible(getOption("zellkonverter.verbose"))
}

.get_verbose <- function(envir) {

    verbose <- envir$verbose

    if (is.null(verbose)) {
        verbose <- getOption("zellkonverter.verbose")
    }

    isTRUE(verbose)
}

.ui_rule <- function(msg, ...) {

    envir <- parent.frame()

    if (.get_verbose(envir)) {
        cli::cli_rule(msg, ..., .envir = envir)
    }
}

.ui_info <- function(msg, ...) {

    envir <- parent.frame()

    if (.get_verbose(envir)) {
        cli::cli_alert_info(msg, ..., .envir = envir)
    }
}

.ui_warn <- function(msg, warn = TRUE, ...) {

    envir <- parent.frame()

    msg <- cli::format_message(msg, .envir = envir)

    if (.get_verbose(envir)) {
        cli::cli_alert_warning(msg, ..., .envir = envir)
    }

    if (warn) {
        warning(msg, call. = FALSE)
    }
}

.ui_step <- function(msg, ...) {

    envir <- parent.frame()

    if (.get_verbose(envir)) {
        cli::cli_progress_step(msg, ..., .envir = envir)
    }
}

.ui_process <- function(msg, ...) {

    envir <- parent.frame()

    if (.get_verbose(envir)) {
        cli::cli_process_start(msg, ..., .envir = envir)
    }
}

.ui_process_done <- function(...) {

    envir <- parent.frame()

    if (.get_verbose(envir)) {
        cli::cli_process_done(..., .envir = envir)
    }
}

.trim_path <- function(path, n = 40) {

    path_split <- .split_path(path)

    for (level in seq_along(path_split)) {
        trimmed_path <- do.call(file.path, as.list(path_split))
        trimmed_path <- gsub("^//", "/", trimmed_path)
        if (nchar(trimmed_path) <= n) {
            break
        } else if (nchar(path_split[level]) >= 3) {
            path_split[level] <- "..."
        }
    }

    return(trimmed_path)
}

.split_path <- function(path) {
    if (dirname(path) != path) {
        path <- c(.split_path(dirname(path)), basename(path))
    }

    return(path)
}
