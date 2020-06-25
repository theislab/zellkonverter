# zellkonverter (development version)

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
