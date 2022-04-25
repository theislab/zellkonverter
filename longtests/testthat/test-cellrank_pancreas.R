library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
file <- bfcrpath(cache, "https://figshare.com/ndownloader/files/30683438")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X", "Ms", "Mu", "fit_t", "fit_tau", "fit_tau_", "spliced",
               "unspliced", "velocity", "velocity_u"),
    colData = c("day", "proliferation", "G2M_score", "S_score", "phase",
                "clusters_coarse", "clusters", "clusters_fine", "louvain_Alpha",
                "louvain_Beta", "palantir_pseudotime", "initial_size_spliced",
                "initial_size_unspliced", "initial_size", "n_counts",
                "velocity_self_transition", "terminal_states",
                "terminal_states_probs", "initial_states",
                "initial_states_probs", "velocity_pseudotime", "latent_time",
                "dpt_pseudotime"),
    rowData = c("highly_variable_genes", "gene_count_corr", "means",
                "dispersions", "dispersions_norm", "highly_variable", "fit_r2",
                "fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling",
                "fit_std_u", "fit_std_s", "fit_likelihood", "fit_u0", "fit_s0",
                "fit_pval_steady", "fit_steady_u", "fit_steady_s",
                "fit_variance", "fit_alignment_scaling", "velocity_genes",
                "to.Epsilon.corr", "to.Alpha.corr", "to.Beta.corr",
                "to.Epsilon.qval", "to.Alpha.qval", "to.Beta.qval"),
    metadata = c("T_bwd_params", "clusters_colors", "clusters_fine_colors",
                 "clusters_sizes", "diffmap_evals", "eig_bwd", "eig_fwd",
                 "initial_states_colors", "initial_states_names", "iroot",
                 "louvain_Alpha_colors", "louvain_Beta_colors", "neighbors",
                 "paga", "pca", "recover_dynamics", "terminal_states_colors",
                 "terminal_states_names", "to_terminal_states_colors",
                 "to_terminal_states_names", "velocity_graph",
                 "velocity_graph_neg", "velocity_params"),
    redDim = c("X_diffmap", "X_pca", "X_umap", "macrostates_bwd",
               "macrostates_fwd", "to_terminal_states", "velocity_umap"),
    varm = c("PCs", "loss"),
    colPairs = c("T_bwd", "T_fwd", "connectivities", "distances")
)

missing <- list()

test_that("Reading H5AD works", {
    expect_warning(
        {sce <- readH5AD(file)},
        "The names of these selected var columns have been modified"
    )
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
