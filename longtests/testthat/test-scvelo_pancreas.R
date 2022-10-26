library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30595479")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X", "Ms", "Mu", "fit_t", "fit_tau", "fit_tau_", "spliced",
               "unspliced", "variance_velocity", "velocity", "velocity_u"),
    colData = c("clusters_coarse", "clusters", "S_score", "G2M_score",
                "initial_size_spliced", "initial_size_unspliced",
                "initial_size", "n_counts", "velocity_self_transition", "phase",
                "velocity_length", "velocity_confidence",
                "velocity_confidence_transition", "root_cells", "end_points",
                "velocity_pseudotime", "latent_time"),
    rowData = c("highly_variable_genes", "gene_count_corr", "means",
                "dispersions", "dispersions_norm", "highly_variable",
                "velocity_gamma", "velocity_qreg_ratio", "velocity_r2",
                "velocity_genes", "spearmans_score", "velocity_score",
                "fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling",
                "fit_std_u", "fit_std_s", "fit_likelihood", "fit_u0", "fit_s0",
                "fit_pval_steady", "fit_steady_u", "fit_steady_s",
                "fit_variance", "fit_alignment_scaling", "fit_r2"),
    metadata = c("clusters_coarse_colors", "clusters_colors", "clusters_sizes",
                 "day_colors", "neighbors", "paga", "pca",
                 "rank_dynamical_genes", "rank_velocity_genes",
                 "recover_dynamics", "velocity_graph", "velocity_graph_neg",
                 "velocity_params"),
    redDim = c("X_pca", "X_umap", "velocity_umap"),
    varm = c("loss"),
    colPairs = c("connectivities", "distances")
)

missing <- list()

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")
})

sce <- suppressWarnings(readH5AD(file))

test_that("SCE is valid", {
    validateH5ADSCE(sce, names, missing)
})

test_that("Writing H5AD works", {
    writeH5AD(sce, outfile)
    expect_true(file.exists(outfile))
})

test_that("Round trip is as expected", {
    out <- readH5AD(outfile)
    expectSCE(out, sce)
})
