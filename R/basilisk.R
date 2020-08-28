#' AnnData dependencies
#'
#' Vector defining a set of Python dependencies and versions required to operate
#' with AnnData and H5AD files.
#'
#' @details
#' This variable is exposed for use by other package developers who want an easy
#' way to define the dependencies required for creating a Python environment to
#' work with AnnData objects, most typically within a **basilisk** context. For
#' example, we can simply combine this vector with additional dependencies to
#' create a **basilisk** environment with Python package versions that are
#' consistent with those in **zellkonverter**.
#'
#' @format
#' A character vector containing the pinned versions of all Python packages on
#' which AnnData depends.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @examples
#' .AnnDataDependencies
#'
#' @export
#' @name AnnDataDependencies
.AnnDataDependencies <- c(
    "anndata==0.7.4",
    "h5py==2.10.0",
    "hdf5==1.10.5",
    "natsort==7.0.1",
    "numpy==1.18.5",
    "packaging==20.4",
    "pandas==1.0.4",
    "scipy==1.4.1",
    "sqlite==3.30.1"
)

anndata_env <- basilisk::BasiliskEnvironment(
    envname="anndata_env",
    pkgname="zellkonverter",
    packages=.AnnDataDependencies
)
