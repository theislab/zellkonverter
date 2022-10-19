# yoinked from reticulate ->
# https://github.com/rstudio/reticulate/blob/fe0eda154a80b22c0d45e043b74390b73ab8b64e/R/utils.R#L49
yoink <- function(package, symbol) {
    do.call(":::", list(package, symbol))
}
disable_conversion_scope <- yoink("reticulate", "disable_conversion_scope")
