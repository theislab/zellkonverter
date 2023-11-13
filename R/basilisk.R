#' AnnData environment
#'
#' The Python environment used by **zellkonverter** for interfacing with the
#' **anndata** Python library (and H5AD files) is described by the dependencies
#' in returned by `AnnDataDependencies()`. The `zellkonverterAnnDataEnv()`
#' functions returns the [basilisk::BasiliskEnvironment()] containing these
#' dependencies used by **zellkonverter**. Allowed versions of **anndata** are
#' available in `.AnnDataVersions`.
#'
#' @details
#'
#' ## Using Python environments
#'
#' When a **zellkonverter** is first run a conda environment containing all of
#' the necessary dependencies for that version with be instantiated. This will
#' not be performed on any subsequent run or if any other **zellkonverter**
#' function has been run prior with the same environment version.
#'
#' By default the **zellkonverter** conda environment will become the shared R
#' Python environment if one does not already exist. When one does exist (for
#' example when a **zellkonverter** function has already been run using a
#' a different environment version) then a separate environment will be used.
#' See [basilisk::setBasiliskShared()] for more information on this behaviour.
#' Note the when the environment is not shared progress messages are lost.
#'
#' ## Development
#'
#' The `AnnDataDependencies()` function is exposed for use by other package
#' developers who want an easy way to define the dependencies required for
#' creating a Python environment to work with AnnData objects, most typically
#' within a **basilisk** context. For example, we can simply combine this
#' vector with additional dependencies to create a **basilisk** environment with
#' Python package versions that are consistent with those in **zellkonverter**.
#'
#' If you want to run code in the exact environment used by **zellkonverter**
#' this can be done using `zellkonverterAnnDataEnv()` in combination with
#' [basilisk::basiliskStart()] and/or [basilisk::basiliskRun()]. Please refer to
#' the **basilisk** documentation for more information on using these
#' environments.
#'
#' @author Luke Zappia
#' @author Aaron Lun
#'
#' @examples
#' .AnnDataVersions
#'
#' AnnDataDependencies()
#' AnnDataDependencies(version = "0.7.6")
#'
#' cl <- basilisk::basiliskStart(zellkonverterAnnDataEnv())
#' anndata <- reticulate::import("anndata")
#' basilisk::basiliskStop(cl)
#' @name AnnData-Environment
#' @rdname AnnData-Environment
NULL

#' @rdname AnnData-Environment
#'
#' @format
#' For `.AnnDataVersions` a character vector containing allowed **anndata**
#' version strings.
#'
#' @export
.AnnDataVersions <- c("0.10.2", "0.9.2", "0.8.0", "0.7.6")

#' @rdname AnnData-Environment
#'
#' @param version A string giving the version of the **anndata** Python library
#' to use. Allowed values are available in `.AnnDataVersions`. By default the
#' latest version is used.
#'
#' @returns
#' For `AnnDataDependencies` a character vector containing the pinned versions
#' of all Python packages to be used by `zellkonverterAnnDataEnv()`.
#'
#' @export
AnnDataDependencies <- function(version = .AnnDataVersions) {

    version <- match.arg(version)

    switch (
        version,
        "0.7.6" = c(
            "anndata==0.7.6",
            "h5py==3.2.1",
            "hdf5==1.10.6",
            "natsort==7.1.1",
            "numpy==1.20.2",
            "packaging==20.9",
            "pandas==1.2.4",
            "python==3.7.10",
            "scipy==1.6.3",
            "sqlite==3.35.5"
        ),
        "0.8.0" = c(
            "anndata==0.8.0",
            "h5py==3.6.0",
            "hdf5==1.12.1",
            "natsort==8.1.0",
            "numpy==1.22.3",
            "packaging==21.3",
            "pandas==1.4.2",
            "python==3.8.13",
            "scipy==1.7.3",
            "sqlite==3.38.2"
        ),
        "0.9.2" = c(
            "anndata==0.9.2",
            "h5py==3.9.0",
            "hdf5==1.14.2",
            "natsort==8.4.0",
            "numpy==1.26.0",
            "packaging==23.2",
            "pandas==2.1.1",
            "python==3.11.5",
            "scipy==1.11.3"
        ),
        "0.10.2" = c(
            "anndata==0.10.2",
            "h5py==3.10.0",
            "hdf5==1.14.2",
            "natsort==8.4.0",
            "numpy==1.26.0",
            "packaging==23.2",
            "pandas==2.1.1",
            "python==3.11.5",
            "scipy==1.11.3"
        )
    )
}

#' @rdname AnnData-Environment
#'
#' @return
#' For `zellkonverterAnnDataEnv` a [basilisk::BasiliskEnvironment()] containing
#'  **zellkonverter**'s AnnData Python environment.
#'
#' @include ui.R
#' @export
zellkonverterAnnDataEnv <- function(version = .AnnDataVersions) {

    version <- match.arg(version)

    basilisk::BasiliskEnvironment(
        envname = paste0("zellkonverterAnnDataEnv-", version),
        pkgname = "zellkonverter",
        packages = AnnDataDependencies(version)
    )
}

# Instantiate environments so they can be found by
# `basilisk::configureBasiliskEnv()` when `BASILISK_USE_SYSTEM_DIR=1`.
# See https://github.com/theislab/zellkonverter/issues/66.
anndata_env_0.7.6 <- zellkonverterAnnDataEnv(version = "0.7.6")
anndata_env_0.8.0 <- zellkonverterAnnDataEnv(version = "0.8.0")
anndata_env_0.9.2 <- zellkonverterAnnDataEnv(version = "0.9.2")
anndata_env_0.10.2 <- zellkonverterAnnDataEnv(version = "0.10.2")
