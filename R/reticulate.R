#' Convert between Python and R objects
#'
#' @param x A Python object.
#'
#' @return An \R object, as converted from the Python object.
#'
#' @details
#' These functions are extensions of the default conversion functions in the
#' `reticulate` package for the following reasons:
#'
#' - `numpy.ndarray` - Handle conversion of **numpy** recarrays
#'
#' @author Luke Zappia
#'
#' @seealso
#' [reticulate::py_to_r()] for the base `reticulate` functions
#'
#' @name r-py-conversion
#' @export
py_to_r.numpy.ndarray <- function(x) {
    disable_conversion_scope(x)

    # Suggested method to detect recarrays from
    # https://stackoverflow.com/a/62491135/4384120
    if (!is.null(py_to_r(x$dtype$names))) {
        # Convert via pandas DataFrame as suggested here
        # https://stackoverflow.com/a/60614003/4384120
        # Not as efficient but less messing around with types
        pandas <- import("pandas", convert = FALSE)
        out <- tryCatch(
            {
                x <- pandas$DataFrame(x)$to_numpy()
                py_to_r(x)
            }, error = function(err) {
                stop("Failed to convert recarray with error: ", err$message,
                     call. = FALSE)
            }
        )
        return(out)
    }

    # No special handler found, delegate to next method
    NextMethod()
}
