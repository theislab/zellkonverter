# zellkonverter

<!-- badges: start -->
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Bioconductor since](https://bioconductor.org/shields/years-in-bioc/zellkonverter.svg)](https://bioconductor.org/packages/devel/bioc/html/zellkonverter.html)
[![Bioconductor rank](https://bioconductor.org/shields/downloads/devel/zellkonverter.svg)](https://bioconductor.org/packages/devel/bioc/html/zellkonverter.html)
[![Bioconductor update](https://bioconductor.org/shields/lastcommit/devel/bioc/zellkonverter.svg)](https://bioconductor.org/packages/devel/bioc/html/zellkonverter.html)
[![Bioconductor build](https://bioconductor.org/shields/build/devel/bioc/zellkonverter.svg)](https://bioconductor.org/packages/devel/bioc/html/zellkonverter.html)
[![R build status](https://github.com/theislab/zellkonverter/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/theislab/zellkonverter/actions)
[![Codecov test coverage](https://codecov.io/gh/theislab/zellkonverter/branch/master/graph/badge.svg)](https://codecov.io/gh/theislab/zellkonverter?branch=master)
<!-- badges: end -->

**zellkonverter** is a small package for converting between SingleCellExperiment
objects and alternative objects for storing single-cell RNA-sequencing data
(such as AnnData). It is built on top of the [**basilisk**][basilisk] package.

For documentation see please refer to [Bioconductor][bioc]. Development
documentation is also available on [Bioconductor devel][bioc-devel] or the
[pkgdown site][pkgdown].

## Installation

**zellkonverter** can be installed from Bioconductor using the **BiocManager**
package:

```r
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("zellkonverter")
```

## Code of Conduct
  
Please note that the zellkonverter project is released with a
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

[basilisk]: https://www.bioconductor.org/packages/basilisk/ "basilisk on Bioconductor"
[bioc]: https://bioconductor.org/packages/zellkonverter/ "zellkonverter on Bioconductor"
[bioc-devel]: https://bioconductor.org/packages/devel/bioc/html/zellkonverter.html "zellkonverter on Bioconductor devel"
[pkgdown]: https://theislab.github.io/zellkonverter/ "zellkonverter pkgdown site"
