\name{NEWS}
\title{News for Package \pkg{zellkonverter}}

\section{Version 1.12.0, Bioconductor 3.18 Release (October 2023)}{
  \subsection{Major changes}{
  \itemize{
    \item{
      Add environments for \bold{anndata} v0.9.2 and v0.10.2. Version 0.10.20 is
      now the default envrionment for the Python reader/writer.
    }
  }}
  \subsection{Minor changes}{
  \itemize{
    \item{
      Changes for compatibility with \bold{rhdf5} v2.45.1 including enum types
      that simplifies reading of nullable types in the native R reader
    }
    \item{
      Dimensions are now passed correctly when converting the \code{raw} slot
    }
    \item{
      Backed sparse matrices are now converted in \code{AnnData2SCE()}
    }
  }}
}

\section{Version 1.10.0, Bioconductor 3.17 Release (April 2023)}{
  \subsection{Major changes}{
  \itemize{
    \item{
      Add compatibility with the \bold{anndata} v0.8 H5AD format to the the
      native R writer (By @jackkamm and @mtmorgan)
    }
    \item{
      Add functions for converting \bold{pandas} arrays used by \bold{anndata}
      when arrays have missing values
    }
  }}
  \subsection{Minor changes}{
  \itemize{
    \item{
      Add Robrecht Cannoodt and Jack Kamm as contributors!
    }
    \item{
      Minor adjustments to tests to match reader changes
    }
  }}
}

\section{Version 1.8.0, Bioconductor 3.16 Release (October 2022)}{
  \subsection{Major changes}{
  \itemize{
    \item{
      Improve compatibility with the R \bold{anndata} package. This required
      modifying conversion functions so that Python objects are explicitly
      converted rather than relying on automatic conversion.
    }
    \item{
      Added support for \bold{numpy} recarrays. This solves a long-standing
      issue and allows results from \bold{scanpy}'s \code{rank_genes_groups()}
      function to be read.
    }
  }}
  \subsection{Minor changes}{
  \itemize{
    \item{
      The Python version is now pinned in the \bold{anndata} v0.7.6 environment
      for compatibility with changes in \bold{basilisk}
    }
    \item{
      Instatiate Python environments so they can be properly picked up by
      \code{basilisk::configureBasiliskEnv()}
    }
    \item{
      Allow missing obs/var names when \code{use_hdf5 = TRUE}
    }
    \item{
      Minor changes to the UI functions for compatibility with \bold{cli} v3.4.0
    }
    \item{
      Minor changes for compatibility with \bold{Matrix} v1.4-2
    }
    \item{
      Improvements to the UI for warnings
    }
    \item{
      Updates and improvments to tests
    }
  }}
}

\section{Version 1.6.0, Bioconductor 3.15 Release (April 2022)}{
  \subsection{Major changes}{
  \itemize{
    \item{
      Added support for multiple \bold{basilisk} environments with different
      \bold{anndata} versions. Users can now specify the environment to use with
      options in \code{readH5AD()} and \code{writeH5AD()}. To faciliate this
      some exported objects where converted to functions but this should only
      effect developers.
    }
    \item{
      Updated the default environment to use \bold{anndata} v0.8.0. This is a
      major update and files written with v0.8.0 cannot be read by previous
      \bold{anndata} versions. This was the motivation for supporting multiple
      environments and users can select the previous environment with
      \bold{anndata} v0.7.6 if compatibility is required.
    }
    \item{
      Standardise naming in \code{AnnData2SCE()}. Column names of data frames
      and names of list items will now be modified to match R conventions
      (according to \code{make.names()}). When this happens a warning will be
      issued listing the modifications. This makes sure than everything in the
      created \code{SingleCellExperiment} is accessible.
    }
  }}
  \subsection{Minor changes}{
  \itemize{
    \item{
      Allow \code{data.frame}'s stored in \code{varm} to be converted in
      \code{SCE2AnnData()}
    }
    \item{
      Minor updates to the vignette and other documentation.
    }
    \item{
      Updates to tests to match the changes above.
    }
  }}
}

\section{Version 1.4.0, Bioconductor 3.14 Release (October 2021)}{
  \itemize{
    \item{
      Add arguments to control how slots are converted in \code{AnnData2SCE()}
      and \code{SCE2AnnData()}. Each slot can now be fully converted, skipped
      entirely or only selected items converted.
    }
    \item{
      Add support for converting the \code{raw} slot to an \code{altExp} in
      \code{AnnData2SCE()}
    }
    \item{
      Add recursive conversion of lists in \code{AnnData2SCE()}
    }
    \item{
      Add progress messages to various functions. These can be controlled by
      function arguments or a global variable.
    }
    \item{
      Add long tests for various public datasets. This should help to make the
      package more robust
    }
    \item{
      Fix bug in converting \code{dgRMatrix} sparse matrices
    }
    \item{
      Correctly handle \code{DataFrame} objects stored in \code{adata.obsm}
    }
  }
}

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
