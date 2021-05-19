\name{NEWS}
\title{News for Package \pkg{zellkonverter}}

\section{Version 1.2.0, Bioconductor 3.13 Release (May 2021)}{
  \itemize{
    \item{
      Update \strong{anndata} and other Python dependencies, now using
      \strong{anndata} v0.7.6
    }
    \item{
      Improved conversion checks for all slots in \code{AnnData2SCE()}
    }
    \item{
      Enable return conversion of the \code{varm} slot in \code{AnnData2SCE()}
    }
    \item{
      Avoid converting \code{obsp} and \code{varp} to dense matrices in
      \code{AnnData2SCE()}
    }
    \item{
      \code{AnnData2SCE()} should now always return \code{dgCMatrix} matrices
      when assays are sparse
    }
    \item{
      More consistent conversion of \code{metadata} to \code{uns} in
      \code{SCE2AnnData()}
    }
    \item{
      Handle conversion of list columns in \code{colData} and \code{rowData} in
      \code{SCE2AnnData()}
    }
    \item{
      Better support for converting \strong{anndata} \code{SparseDataset} arrays
    }
    \item{
      Improved support for conversion of HDF5 backed \code{AnnData} objects
    }
    \item{
      Better support for writing \code{DelayedArray} assays in
      \code{writeH5AD()}
    }
    \item{
      Store \code{X_name} in \code{AnnData2SCE()} for use by
      \code{SCE2AnnData()} and add an \code{X_name} argument to
      \code{AnnData2SCE()} and \code{readH5AD()}
    }
    \item{
      Add a \code{compression} argument to \code{writeH5AD()}
    }
    \item{
      Add an experimental native R reader to \code{readH5AD()}
    }
    \item{
      Export \code{zellkonverterAnnDataEnv} for use by other packages
    }
  }
}

\section{Version 1.0.0, Bioconductor 3.12 Release (October 2020)}{
  \itemize{
    \item{Accepted into Bioconductor for Release 3.12}
    \item{
      zellkonverter provides methods to convert between Python AnnData objects
      and SingleCellExperiment objects. These are primarily intended for use by
      downstream Bioconductor packages that wrap Python methods for single-cell
      data analysis. It also includes functions to read and write H5AD files
      used for saving AnnData objects to disk.
    }
  }
}
