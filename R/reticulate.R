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
#' - `pandas.core.arrays.masked.BaseMaskedArray` - Handle conversion of
#'   **pandas** arrays (used when by `AnnData` objects when there are missing
#'   values)
#' - `pandas.core.arrays.categorical.Categorical` - Handle conversion of
#'   **pandas** categorical arrays
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

#' @export
py_to_r.pandas.core.arrays.masked.BaseMaskedArray <- function(x) {
    disable_conversion_scope(x)

    if (is(x, "pandas.core.arrays.boolean.BooleanArray")) {
        dtype <- "bool"
        fill <- FALSE
    } else if (is(x, "pandas.core.arrays.integer.IntegerArray")) {
        dtype <- "int"
        fill <- 0L
    } else if (is(x, "pandas.core.arrays.floating.FloatingArray")) {
        dtype <- "float"
        fill <- 0.0
    } else if (is(x , "pandas.core.arrays.string_.StringArray")) {
        dtype <- "str"
        fill <- ""
    } else {
        stop(
            "No conversion exists for this Pandas array type: ",
            paste(class(x), collapse = ", ")
        )
    }

    # Record which values should be NA
    is_na <- reticulate::py_to_r(x$isna())

    # Fill NA values with a dummy
    x <- x$fillna(value = fill)

    # Convert to numpy array and then to R using default conversion
    x <- x$to_numpy()$astype(dtype)
    x <- reticulate::py_to_r(x)

    # Restore the NA values
    x[is_na] <- NA

    return(x)
}

#' @export
py_to_r.pandas.core.arrays.categorical.Categorical <- function(x) {
    disable_conversion_scope(x)

    # Get the category levels
    cats <- reticulate::py_to_r(x$categories$to_list())

    # Record which values should be NA
    is_na <- reticulate::py_to_r(x$isna())

    # Fill NA values with a dummy
    x <- x$fillna(value = cats[1])

    # Convert to list and then to R using default conversion
    x <- x$tolist()
    x <- reticulate::py_to_r(x)

    # Restore the NA values
    x[is_na] <- NA

    # Convert to factor
    x <- factor(x, levels = cats)

    return(x)
}
