# This script was used to create the `example_anndata.h5ad` file.
# This file contains an example AnnData object for use in examples and tests.
# A dataset with 200 cells and 500 genes was generated using the Splat simulation in the Splatter package.
# A Python AnnData object was created using this data (via reticulate) and run through a standard Scanpy analysis workflow to populate the various slots.
# The file object was then saved to disk as a .h5ad file.
#
# Key package versions:
#
# splatter   v1.14.0
# reticulate v1.18
# scanpy     v1.5.1
# anndata    v0.7.4

library(splatter)
library(reticulate)

mini_sim <- splatSimulateGroups(batchCells = 200, nGenes = 500, lib.loc = 8,
                                group.prob = c(0.5, 0.5), seed = 1)

anndata <- import("anndata")
scanpy  <- import("scanpy")

adata <- anndata$AnnData(t(counts(mini_sim)))
adata$obs_names <- colnames(mini_sim)
adata$var_names <- rownames(mini_sim)
adata$layers <- list(counts = t(counts(mini_sim)))

scanpy$pp$filter_genes(adata, min_counts = 10)
scanpy$pp$normalize_total(adata, target_sum = 1e4)
scanpy$pp$log1p(adata)
scanpy$pp$highly_variable_genes(adata)
scanpy$tl$pca(adata, svd_solver = "arpack")
scanpy$pp$neighbors(adata, n_pcs = 10L)
scanpy$tl$umap(adata)
scanpy$tl$louvain(adata)
scanpy$tl$rank_genes_groups(adata, "louvain")

adata$write_h5ad("example_anndata.h5ad")
