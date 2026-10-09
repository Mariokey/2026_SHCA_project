# 2026_SHCA_project

Analysis pipelines for a spatially resolved multi-organ single-cell atlas and conserved
CNTN4+ neuroepithelial cell states.

---

## What this repository is

This is an **analysis record**, not a turnkey runnable package. Each directory
corresponds to one section of the Methods, and contains the code that produced
the corresponding results and figures.

The scripts read intermediate objects (Seurat `.rds`, AnnData `.h5ad`) from the
original analysis workspace. Those objects are **not distributed here** — the
atlas comprises 115 samples and is far too large for a code repository. Paths in
the scripts have therefore been normalised to a small set of conventional
locations (see [Data layout](#data-layout)) rather than the absolute paths of
the machine the analysis ran on.

To reuse a script you will need to point it at your own objects, or ask the
authors for access to the intermediates.

## Repository layout

| Directory | Methods section | Contents |
|---|---|---|
| `01_preprocessing_annotation/` | *Preprocessing of single-cell RNA-sequencing data*, *Clustering and cell type annotation* | scanpy QC/merge/annotate pipeline (`preprocess_pipeline.py`); Seurat label curation per compartment |
| `02_gene_score/` | *Gene score analysis* | `AddModuleScore` signature scoring; CNTN4 vs neuron-receptor score correlation |
| `03_deg_function_enrichment/` | *Function annotation* | DEG detection per cluster and per location x subtype; GO/KEGG enrichment |
| `04_cell_modules_covarnet/` | *Identification of cellular modules* | CoVarNet cellular co-occurrence modules; hdWGCNA gene modules |
| `05_nmf/` | *NMF analysis* | NMF factorisation, rank survey, meta-program calling and scoring |
| `06_celltype_proportion_roe/` | *Enrichment analysis of cell type proportion* | Ro/e computation and entropy visualisation |
| `07_trajectory/` | *Trajectory analysis* | Monocle2, Slingshot, tradeSeq, macrophage trajectory scoring |
| `08_cellchat_cci/` | *Cell-cell interactions* | CellChat ligand-receptor inference, centrality, network plots |
| `09_dendrogram_comparison/` | *Dendrogram comparison* | Robinson-Foulds distance to the reference organ-developmental tree |
| `10_public_data_scibet/` | *Public scRNA-seq data processing* | scibet projection of 57 public datasets; CNTN4+ epithelial recovery |
| `11_disease_virus_enrichment/` | *Association analysis of disease and virus infection related gene sets* | DisGeNET / PheWAS risk-gene overlap, hypergeometric enrichment |
| `12_figure_assembly/` | *Statistics and reproducibility* | Composite overview figures and curated formal panels |
| `utils/` | — | Shared R helpers and the Python `sc_utils` module |

Within each directory, files are numbered in the order they were run.

## Data layout

The scripts expect the following, relative to the working directory:

```
./data/      intermediate Seurat/AnnData objects produced by earlier steps
./ref/       reference gene lists (protein-coding genes, TF lists),
             marker signature tables, DisGeNET / PheWAS risk gene sets
DATA_DIR     placeholder for large inputs that lived outside the project
             directory; substitute your own path
```

Objects referenced by name in the scripts include the per-compartment annotated
Seurat objects (`01.epi_seuratAnno.*.rds`, `01.immu_seuratAnno.*.rds`,
`01.strom_seuratAnno.*.rds`, `01.endo_seuratAnno.*.rds`,
`01.germ_seuratAnno.*.rds`, `01.all_seuratAnno.*.rds`), the merged atlas object.

## Software environment

**R** (analysis was performed with R 4.3.0; the repository parses cleanly under
R 4.4.3). Principal packages, with versions where the Methods specify one:

- Seurat 4.3.1, Harmony 1.0, clusterProfiler 3.14.3, ComplexHeatmap 2.10.0,
  pheatmap 1.0.12, NMF, hdWGCNA, CoVarNet, slingshot 2.2.1, tradeSeq 1.8.0,
  monocle 2.14.0, CellChat 1.6.1, scibet, Treedist 2.9.1, HGNChelper 0.8.1,
  AnnotationDbi 1.148.0, org.Hs.eg.db 3.10.0, ggplot2 3.3.6
- The original sessions sourced a personal library-loading helper. That file is
  not included; load the packages above directly, or adapt `utils/00_env_setup.R`
  to your own environment.
- A small number of helper functions (`modifyAno`, `runFig`, `subsetID2`,
  `MultiMrun`, `scPalette2`, `rowPercents`/`colPercents`, `DimPlot_idx`) come
  from the authors' internal R utility library and are called but not defined
  here. They are thin wrappers; their behaviour is documented inline where used.

**Python** (scanpy-based preprocessing):

- scanpy >= 1.9, numpy, pandas, seaborn, matplotlib, anndata
- scrublet is used when installed; otherwise the pipeline falls back to
  scanpy's built-in `sc.external.pp.scrublet`
- harmonypy (only when Harmony batch correction is enabled)

## Running the preprocessing pipeline

`01_preprocessing_annotation/preprocess_pipeline.py` is the one entry point that
is written to run end to end on your own data:

```bash
python 01_preprocessing_annotation/preprocess_pipeline.py \
    --h5ad-dir ./data/per_sample_h5ad \
    --meta     ./data/sample_sheet.csv \
    --markers  ./ref/marker_signatures.tsv \
    --out      ./data/merged_annotated.h5ad
```

Inputs:

- `--h5ad-dir` — one `.h5ad` per sample holding **raw counts**
- `--meta` — sample sheet with a `SampleID` column; extra columns are joined
  into `obs`
- `--markers` — long-format TSV with `celltype` and `gene` columns

## Reproducing figures

As stated in the Methods, all figures can be reproduced with the Seurat
functions `DimPlot`, `FeaturePlot`, `DotPlot` and `VlnPlot`, with heatmaps built
in `ComplexHeatmap`. `12_figure_assembly/` collects the curated panel code;
`utils/helper_plotting.R` holds the shared plotting wrappers those panels call.

