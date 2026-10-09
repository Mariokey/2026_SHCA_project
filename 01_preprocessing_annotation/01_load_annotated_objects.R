# ============================================================================
# Loading the per-compartment annotated Seurat objects
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Preprocessing of single-cell RNA-sequencing data"
# Source     : extracted from Analysis240525.r (L85-106)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
###20250227
seu_all91w=h5ad2seu("../data/AllCombined_Filter.acdata_qc_h5.91w_harmony0418.h5ad", dataname = "01.all91wlarge", savedir = "./")
print(colnames(seu_all91w@meta.data))

major_order=c(
  "Epithelial Cells",
  "Endothelial Cells",
  "Stromal Cells",
  "Immune Cells",
  "Germ Cells",
  "Prolif.")

seu_immu=readRDS("./01.immu_seuratAnno.v0525.rds")
seu_epi=readRDS("./01.epi_seuratAnno.v0525.rds")
seu_strom=readRDS("./01.strom_seuratAnno.v0525.rds")
seu_endo=readRDS("./01.endo_seuratAnno.v0525.rds")
seu_germ=readRDS("./01.germ_seuratAnno.v0525.rds")
seu_prof=readRDS("./01.prof_seuratAnno.v0525.rds")
seu_all=readRDS("./01.all_seuratAnno.v0525.rds")
