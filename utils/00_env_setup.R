# ============================================================================
# Environment setup, plotting themes and output directories
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Source     : extracted from Analysis240525.r (L1-18)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
# NOTE: the original session sourced a personal package-loading helper here
#       (library(...) calls for Seurat, ggplot2, ComplexHeatmap, pheatmap, etc.).
#       Load the packages listed in the repository README instead.
setwd("./data/")
library(reticulate)
reticulate::use_condaenv("SCS")

blank= theme(axis.title=element_blank(),
             axis.text=element_blank(),
             axis.ticks=element_blank())
myfont= theme(text=element_text(size=14,  family="Helvetica", colour ="black"))
NoTitle=theme(plot.title = element_blank())
ThemeBox=theme_bw()+theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

OutputDir="./AnaRes20240524/"
