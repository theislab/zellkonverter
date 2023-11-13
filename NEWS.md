# zellkonverter 1.12.0

* Bioconductor 3.18, October 2023

## zellkonverter 1.12.1 (2023-11-13)

* Fix the **anndata** v0.10.2 environment instatiation (Fixes #103)
* Fix a type in the AnnData Conversion docs (Fixes #100)

## zellkonverter 1.11.4 (2023-10-16)

* Add environment for **anndata** v0.10.2

## zellkonverter 1.11.3 (2023-10-2)

* Add environment for **anndata** v0.9.2

## zellkonverter 1.11.2 (2023-08-28)

* Changes for compatibility with **{rhdf5}** v2.45.1
  * Support for enum types that simplifies reading of nullable types in the
    native R reader

## zellkonverter 1.11.1 (2023-05-23)

* Pass correct dimensions when converting `raw` (Fixes #96)
* Convert **anndata** backed sparse matrices in `AnnData2SCE()` (Fixes #96)

## zellkonverter 1.11.0 (2023-04-26)

* Bioconductor 3.18 devel

# zellkonverter 1.10.0 (2023-04-26)

* Bioconductor 3.17, April 2023

## zellkonverter 1.10.1 (2023-05-23)

* Pass correct dimensions when converting `raw` (Fixes #96)
* Convert **anndata** backed sparse matrices in `AnnData2SCE()` (Fixes #96)

## zellkonverter 1.9.3 (2023-04-06)

* Add functions for converting **pandas** arrays used by **anndata** when
  arrays have missing values (Fixes #87)
* Read the correct index names in the R reader (PR #93 mtmorgan)
* Adjust tests to match reader changes

## zellkonverter 1.9.2 (2023-03-28)

* Add @rcannood as a contributor (PR #90 @rcannood, fixes #88)

## zellkonverter 1.9.1 (2023-03-14)

* Add compatibility with the **anndata** v0.8 H5AD format to the the native R
  writer (PR #86 @jackkamm, fixes #78)

## zellkonverter 1.9.0 (2022-11-02)

* Bioconductor 3.17 devel

# zellkonverter 1.8.0 (2022-11-02)

* Bioconductor 3.16, November 2022

## zellkonverter 1.7.8 (2022-10-04)

* Improve compatibility with the R **{anndata}** package (PR #76 @rcannood,
  fixes #75)
  * Python objects are now explicitly converted rather than relying on automatic
    conversion
  * Other minor modifications for compatibility
* Added support for **numpy** recarrays (dtype number 20) (PR #81, fixes #45,
  #28)
  * Added a new `py_to_r.numpy.ndarray()` function which extends the default
    **{reticulate}** function
* Improvements to warnings
* Improvements and updates to tests

## zellkonverter 1.7.7 (2022-10-04)

* Pin **python** version to 3.7.10 in **anndata** v0.7.6 environment (3.7.12
  was not compatible with other dependencies)

## zellkonverter 1.7.6 (2022-09-29)

* Pin **python** version to 3.7.12 in **anndata** v0.7.6 environment to match
  **{basilisk}** changes

## zellkonverter 1.7.5 (2022-09-13)

* Minor changes for compatibility with **{cli}** v3.4.0
  * Added tests for `verbose=TRUE` 

## zellkonverter 1.7.4 (2022-08-17)

* Minor changes for compatibility with the upcoming **{Matrix}** 1.4-2 release

## zellkonverter 1.7.3 (2022-06-23)

* Move verbose from `zellkonverterAnnDataEnv()` (Fixes #66)

## zellkonverter 1.7.2 (2022-06-09)

* Instantiate environments for `basilisk::configureBasiliskEnv()` (Fixes #66)
* Allow missing obs/var names when `use_hdf5 = TRUE` (Fixes #65)

## zellkonverter 1.7.1 (2022-05-17)

* Fix bug in long tests

## zellkonverter 1.7.0 (2022-04-27)

* Bioconductor 3.16 devel

# zellkonverter 1.6.0 (2022-04-27)

* Bioconductor 3.15, April 2022

## zellkonverter 1.6.5 (2022-09-13)

* Minor changes for compatibility with **{cli}** v3.4.0
  * Added tests for `verbose=TRUE` 

## zellkonverter 1.6.4 (2022-08-17)

* Minor changes for compatibility with the upcoming **{Matrix}** 1.4-2 release

## zellkonverter 1.6.3 (2022-06-23)

* Move verbose from `zellkonverterAnnDataEnv()` (Fixes #66)

## zellkonverter 1.6.2 (2022-06-09)

* Instantiate environments for `basilisk::configureBasiliskEnv()` (Fixes #66)
* Allow missing obs/var names when `use_hdf5 = TRUE` (Fixes #65)

## zellkonverter 1.6.1 (2022-05-17)

* Fix bug in long tests

## zellkonverter 1.5.4 (2022-04-25)

* Fix progress messages in `.convert_anndata_df()`
* Allow `data.frames` in `varm` in `SCE2AnnData()`
* Standardise `uns` names to match R conventions in `AnnData2SCE()`
* Adjust long tests

## zellkonverter 1.5.3 (2022-04-19)

* Reduce **scipy** version to 1.7.3
  * **scipy** >= 1.8.0 is incompatible with **{reticulate}** <= 1.24 (see
    https://github.com/rstudio/reticulate/pull/1173)
* Add GTEX 8 tissues dataset to long tests (see #58)

## zellkonverter 1.5.2 (2022-04-17)

* Update the default Python environment to use **anndata** v0.8.0
    * **anndata** 0.8.0
    * **h5py** 3.6.0
    * **hdf5** 1.12.1
    * **natsort** 8.1.0
    * **numpy** 1.22.3
    * **packaging** 21.3
    * **pandas** 1.4.2
    * **python** 3.8.13
    * **scipy** 1.8.0
    * **sqlite** 3.38.2
* Add options to choose Python environments with different versions of
  **anndata**
  * To facilitate this `zellkonverterAnnDataEnv()` and `AnnDataDependencies()`
    are new functions rather than variables
  * Added a new `.AnnDataVersions` variable which stores the available
    **anndata** versions
  * Updates to the vignette and function documentation explaining this option

## zellkonverter 1.5.1 (2022-03-21)

* Modify how Pandas DataFrames are converted to R
  * Columns should now use R approved names with a warning when changes are
    made

## zellkonverter 1.5.0 (2021-10-27)

* Bioconductor 3.15 devel

# zellkonverter 1.4.0 (2021-10-27)

* Bioconductor 3.14, October 2021

## zellkonverter 1.3.3 (2021-10-20)

* Add progress messages to various functions
  * Can be controlled by function arguments or a global variable
* Split `konverter.R` into two files (`AnnData2SCE.R` and `SCE2AnnData.R`)
* Add arguments to control how slots are converted in `AnnData2SCE()` and
  `SCE2AnnData()` (Fixes #47)
  * Each slot can now be fully converted, skipped entirely or only selected
    items converted.
* Add support for converting the `raw` slot to an `altExp` in `AnnData2SCE()` 
  (Fixes #53, fixes #57)

## zellkonverter 1.3.2 (2021-09-09)

* Add recursive conversion of lists in `AnnData2SCE()`
* Correctly handle `DataFrame` objects stored in `adata.obsm`
* Remove **pandas** indexes from converted `DataFrame` objects
* Add functions for validating `SingleCellExperiment` objects (for testing)
* Add long tests for various public datasets

## zellkonverter 1.3.1 (2021-06-22)

* Fix bug in converting `dgRMatrix` sparse matrices (Fixes #55)

## zellkonverter 1.3.0 (2021-05-20)

* Bioconductor 3.14 devel

# zellkonverter 1.2.0 (2021-05-20)

* Bioconductor 3.13, May 2021

## zellkonverter 1.2.1 (2021-06-22)

* Fix bug in converting `dgRMatrix` sparse matrices (Fixes #55)

## zellkonverter 1.1.11 (2021-05-19)

* Add experimental native R reader to `readH5AD()`

## zellkonverter 1.1.10 (2021-05-18)

* Update NEWS for release

## zellkonverter 1.1.9 (2021-05-12)

* `AnnData2SCE()` no longer returns `dgRMatrix` sparse matrices (Fixes #34)

## zellkonverter 1.1.8 (2021-05-03)

* Add conversion checks to all slots in `AnnData2SCE()` (See #45)
* Enable return conversion of `varm` in `SCE2AnnData()` (Fixes #43)
* Store `X_name` in `AnnData2SCE()` for use by `SCE2AnnData()` and add an
  `X_name` argument to `AnnData2SCE()` and `readH5AD()` (Fixes #7)

## zellkonverter 1.1.7 (2021-04-30)

* Add `compression` argument to `writeH5AD()` (Fixes #49)
* Update **anndata** Python dependencies, now using **anndata** v0.7.6

## zellkonverter 1.1.6 (2021-04-27)

* Adapt to changes in `HDF5Array::HDF5Array()`

## zellkonverter 1.1.5 (2021-03-05)

* Better support for **anndata** `SparseDataset` arrays (PR #41, Fixes #37,
  Fixes #42)
* More consistent conversion of `metadata` to `uns` in `SCE2AnnData()`
  (Fixes #40)
* Add handling of list columns in `colData` and `rowData` in `SCE2AnnData()`
  (Fixes #26)
* Export `zellkonverterAnnDataEnv` (Fixes #38)

## zellkonverter 1.1.4 (2021-02-18)

* Handle writing **DelayedArray** assays on the R side in `writeH5AD()`
  (PR #35, Fixes #32)

## zellkonverter 1.1.3 (2021-01-22)

* Adjust `SCE2AnnData()` example (Fixes #31)

## zellkonverter 1.1.2 (2020-12-19)

* Improved support for HDF5 backed conversion (PR #29, fixes #13)

## zellkonverter 1.1.1 (2020-12-03)

* Add `example_anndata.h5ad` file to `inst/extdata/` and creation script to `inst/scripts/`
* Improve conversion checks when converting `.uns` to `metadata`
* Avoid converting `obsp` and `varp` to dense matrices

## zellkonverter 1.1.0 (2020-10-28)

* Bioconductor 3.13 devel

# zellkonverter 1.0.0 (2020-10-28)

* Bioconductor 3.12, October 2020

## zellkonverter 1.0.3 (2021-03-08)

* Avoid converting `obsp` and `varp` to dense matrices

## zellkonverter 1.0.2 (2021-01-28)

* Merge remaining commits for HDF5 conversion (fixes #33)

## zellkonverter 1.0.1 (2021-01-26)

* Improved support for HDF5 backed conversion (PR #29, fixes #13, fixes #33)

## zellkonverter 0.99.7 (2020-10-16)

* Update Python dependencies
  * **numpy** 1.18.5 -> 1.19.1
  * **pandas** 1.0.4 -> 1.1.2
  * **scipy** 1.4.1 -> 1.5.2
  * **sqlite** 3.30.1 -> 3.33.0

## zellkonverter 0.99.6 (2020-10-12)

* Document character to factor coercion in `writeH5ad()` (Fixes #6)
* Add `X_name` argument to `writeH5AD()` (Fixes #23)

## zellkonverter 0.99.5 (2020-09-28)

* Tidy NEWS files for Bioconductor release

## zellkonverter 0.99.4 (2020-08-28)

* Bump anndata version to 0.7.4

## zellkonverter 0.99.3 (2020-08-21)

* Document the `krumsiek11.h5ad` file
* Remove the `internal` keyword from the `zellkonverter-package` documentation

## zellkonverter 0.99.2 (2020-08-21)

* Update `.gitignore`

## zellkonverter 0.99.1 (2020-07-15)

* Fix SCE to AnnData map figure in PDF manual
* Use `expect_equal()` instead of `expect_identical()` in `writeH5AD()` sparse
  matrices test
* Edit package title and description

## zellkonverter 0.99.0 (2020-07-10)

* Initial Bioconductor submission

# zellkonverter 0.0.0 (early development version)

## zellkonverter 0.0.0.9017 (2020-07-10)

* Add biocViews to DESCRIPTION
* Edit package description
* Tidy code
* Replace 1:... with `seq_len()`

## zellkonverter 0.0.0.9016 (2020-07-10)

* Add check for **scRNAseq** in examples (Fixes #18)

## zellkonverter 0.0.0.9015 (2020-07-02)

* Skip `AnnData` matrices without a transposable R counterpart
* Only replace skipped matrices when `use_hdf5 = TRUE` in `readH5AD()`
  (Fixes #12)
* Additional tests for sparse matrices

## zellkonverter 0.0.0.9014 (2020-06-30)

* Allow assay skipping when converting from `SingleCellExperiment` to `AnnData`
* Allow skipping of assays that aren't **numpy** friendly in `writeH5AD()`
* Wait for **basilisk** process shutdown to release `.h5ad` file
* Updates to documentation and tests

## zellkonverter 0.0.0.9013 (2020-06-25)

* Improve conversion between `SingleCellExperiment` and `AnnData` (See #8)
  * Convert between `metadata` and `uns` (where objects are compatible)
  * Convert between `rowPairs` and `varp`
  * Convert between `colPairs` and `obsp`
  * Convert from `varm` to `rowData` (but not in reverse)
* Add mapping table to docs

## zellkonverter 0.0.0.9012 (2020-06-19)

* Tidy documentation and code
* Tidy vignette

## zellkonverter 0.0.0.9011 (2020-06-18)

* Support for HDF5Array outputs in `readH5AD()` (Fixes #4)

## zellkonverter 0.0.0.9010 (2020-06-17)

* Avoid checking column names for `colData` and `rowData` in `SCE2AnnData()`
* Make sure that all matrices passes to **{reticulate}** are **numpy** friendly
* Add more tests
* Update vignette front matter

## zellkonverter 0.0.0.9009 (2020-06-15)

* Add vignette

## zellkonverter 0.0.0.9008 (2020-06-12)

* Add examples and improve documentation
* Export `.AnnDataDependencies` for external use

## zellkonverter 0.0.0.9007 (2020-06-11)

* Add `SCE2AnnData()` function
* Add `writeH5AD()` function

## zellkonverter 0.0.0.9006 (2020-06-11)

* Use internal function in `readH5AD()`

## zellkonverter 0.0.0.9005 (2020-06-09)

* Rename `adata2SCE()` to `AnnData2SCE()`
* Remove **{basilisk}** context from `AnnData2SCE()` (See #1)
  * Now uses the calling context

## zellkonverter 0.0.0.9004 (2020-06-09)

* Pin more **AnnData** dependencies (See #1)

## zellkonverter 0.0.0.9003 (2020-06-08)

* Add test `.h5ad` file
* Add test for `readH5AD()`
* Add package man page

## zellkonverter 0.0.0.9002 (2020-06-08)

* Add `adata2SCE()` function
* Add `readH5AD()` function

## zellkonverter 0.0.0.9001 (2020-06-08)

* Add **{basilisk}** infrastructure

## zellkonverter 0.0.0.9000 (2020-06-08)

* Set up package
