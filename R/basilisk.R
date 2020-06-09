anndata_deps <- c(
    "anndata==0.7.3",
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
    packages=anndata_deps
)
