#' AnnData environment
#'
#' The Python environment used by **zellkonverter** for interfacing with the
#' **anndata** Python library (and H5AD files) is described by the dependencies
#' in `.AnnDataDependencies`. The `zellkonverterAnnDataEnv` variable is the
#' [basilisk::BasiliskEnvironment()] used by **zellkonverter**.
#'
#' @details
#' The `.AnnDataDependencies` variable is exposed for use by other package
#' developers who want an easy way to define the dependencies required for
#' creating a Python environment to work with AnnData objects, most typically
#' within a **basilisk** context. For example, we can simply combine this
#' vector with additional dependencies to create a **basilisk** environment with
#' Python package versions that are consistent with those in **zellkonverter**.
#'
#' If you want to run code in the exact environment used by **zellkonverter**
#' this can be done using `zellkonverterAnnDataEnv` in combination with
#' [basilisk::basiliskStart()] and/or [basilisk::basiliskRun()]. Please refer to
#' the **basilisk** documentation for more information on using these
#' environments.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @examples
#' .AnnDataDependencies
#'
#' cl <- basilisk::basiliskStart(zellkonverterAnnDataEnv)
#' anndata <- reticulate::import("anndata")
#' basilisk::basiliskStop(cl)
#'
#' @name AnnData-Environment
#' @rdname AnnData-Environment
NULL

#' @rdname AnnData-Environment
#'
#' @format
#' A character vector containing the pinned versions of all Python packages in
#' `zellkonverterAnnDataEnv`.
#'
#' @export
.AnnDataDependencies <- c(
    "anndata==0.7.4",
    "h5py==2.10.0",
    "hdf5==1.10.5",
    "natsort==7.0.1",
    "numpy==1.19.1",
    "packaging==20.4",
    "pandas==1.1.2",
    "scipy==1.5.2",
    "sqlite==3.33.0"
)

#' @rdname AnnData-Environment
#'
#' @format
#' A [basilisk::BasiliskEnvironment()] containing **zellkonverter**'s AnnData
#' Python environment.
#'
#' @export
zellkonverterAnnDataEnv <- basilisk::BasiliskEnvironment(
    envname  = "zellkonverterAnnDataEnv",
    pkgname  = "zellkonverter",
    packages = .AnnDataDependencies
)
