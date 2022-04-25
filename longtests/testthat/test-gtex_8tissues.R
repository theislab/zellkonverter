library(SingleCellExperiment)
library(BiocFileCache)

cache <- BiocFileCache(ask = FALSE)
# Available from https://www.gtexportal.org/home/datasets
file <- bfcrpath(cache, "https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")
outfile <- tempfile(fileext = ".h5ad")

names <- list(
    assays = c("X", "counts"),
    colData = c(
        "n_genes", "fpr", "tissue", "prep", "individual", "nGenes", "nUMIs",
        "PercentMito", "PercentRibo", "Age_bin", "Sex", "Sample.ID",
        "Participant.ID", "Sample.ID.short",
        "RIN.score.from.PAXgene.tissue.Aliquot",
        "RIN.score.from.Frozen.tissue.Aliquot", "Autolysis.Score",
        "Sample.Ischemic.Time..mins.", "Tissue.Site.Detail", "scrublet",
        "scrublet_score", "barcode", "batch", "n_counts",
        "tissue.individual.prep", "Broad.cell.type", "Granular.cell.type",
        "introns", "junctions", "exons", "sense", "antisense", "intergenic",
        "batch.barcode", "exon_ratio", "intron_ratio", "junction_ratio",
        "log10_nUMIs", "leiden", "leiden_tissue", "Tissue.composition",
        "Cell.types.level.2", "Cell.types.level.3", "Broad.cell.type.numbers",
        "Broad.cell.type..numbers.", "Tissue", "channel"
    ),
    rowData = c(
        "gene_ids", "Chromosome", "Source", "Start", "End", "Strand",
        "gene_name", "gene_source", "gene_biotype", "gene_length",
        "gene_coding_length", "Approved.symbol", "Approved.name", "Status",
        "Previous.symbols", "Alias.symbols", "gene_include", "n_cells"
    ),
    metadata = c(
        "Broad.cell.type..numbers._colors", "Broad.cell.type.numbers_colors",
        "Broad.cell.type_colors", "Broad.cell.type_logregcv_vae_colors",
        "Broad.cell.type_sizes", "Granular.cell.type_colors",
        "Participant.ID_colors", "Sex_colors", "Tissue.composition_colors",
        "Tissue_colors", "dendrogram_..Broad.cell.type..", "leiden",
        "leiden_colors", "leiden_sub_colors", "neighbors", "paga",
        "prep_colors", "tissue_colors", "umap"
    ),
    redDim = c(
        "X_pca", "X_umap", "X_umap_tissue", "X_vae_mean", "X_vae_mean_tissue",
        "X_vae_samples", "X_vae_var"
    ),
    varm = c("spring_leiden_sub"),
    colPairs = c("connectivities", "distances")
)

missing <- list()

test_that("Reading H5AD works", {
    expect_warning(
        {sce <- readH5AD(file)},
        "The names of these selected uns items have been modified"
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
