r_to_py_ifneedbe <- function(x, convert = FALSE) {
  if (inherits(x, "python.builtin.object")) {
    x
  } else {
    r_to_py(x, convert = convert)
  }
}
