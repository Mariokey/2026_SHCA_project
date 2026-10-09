# ============================================================================
# Ro/e cell-type proportion enrichment per compartment and per organ
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Enrichment analysis of cell type proportion"
# Source     : extracted from Analysis240525.r (L546-560, L2735-2745, L7960-7970, L8180-8200)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
#xxx[[1]]$layers[[2]]$show.legend=F
#xxx[[1]]+scale_fill_gradient2(low = "#141464", mid = "white", high = "#78050F", midpoint = -0.05) #, limits=c(-2.5,2.5))
sel_RNA_obj=AverageExpression(Mobj, features = select_major_markers, assays = "RNA", group.by = "subtypes03", return.seurat = TRUE)

pheatmap(scale(t(sel_RNA_obj@assays$RNA$data[select_major_markers,]), center = T, scale = T), name = "Z-score\nof Gene\nExpression",
         color = viridis::viridis(n=100), angle_col = "45",
         cluster_rows = F, cluster_cols = F, scale = F)

Xdata=subset(subset(Mobj, System %in% c("DS")), MajorClass01 %in% c("Immune Cells"))
Xdata$Location_short=droplevels(Xdata$Location_short)
Xdata$subtypes03=droplevels(Xdata$subtypes03)
roveRun(Xdata, Group="Location_short", CT="subtypes03", outdir=OutputDir, mywidth=3.5, nCT=30)

roveRun(Mobj, Group="Location_short", CT="subtypes011", outdir="./", mywidth=20, nCT=107)

###
outdir="./data/R_data/AnaRes20240907/"
dir.create(outdir)
roveRun(inobj, Group="Location_short", CT="MYL_SubClass04",
        outdir=outdir, flip = FALSE, mywidth=8, nCT=15)

xxxn$Group=droplevels(xxxn$Group)
roveRun(xxxn, Group="Group1", CT="ident",         outdir="./rove.pdf", flip = FALSE, mywidth=8, nCT=5)

roveRun(subdataList$strom_fb, Group="Location_short", CT="FB_SubClass02",
        outdir=paste0(outdir,"/FB_"), flip = FALSE, mywidth=8, nCT=15)

xxx2=table(ImmuMData@meta.data[,'SubClass01'],ImmuMData$Location_short)
xxx12=rbind(xxx1,xxx2)
xxx12[xxx12<=10]=0
xxx12=xxx12[rowSums(xxx12)>0,]
xxx12=xxx12[,colSums(xxx12)>0]
yyy=read.table(paste0(dir_img,'LM.rove_main_age.txt'))
yyy=yyy[rownames(yyy) %in% rownames(xxx12) ,]
yyy=yyy[,colnames(yyy) %in% colnames(xxx12)]
plot2=roe_entro_plot(roefile=paste0(dir_img,'LM.rove_main_age.txt'),
                    intable=xxx12, thresh=2)

xxx=FindMarkers(ImmuMData,ident.1 = c("Mph01_AGRP"), group.by = "SubClass01", test.use = "t", only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)

plotF12=roveRun(subdataList$strom_fb, Group="Location_short", CT="FB_SubClass02", outdir=paste0("./FBnew."), flip = FALSE, mywidth=4.5, nCT=13)
plotF2=roe_entro_plot(roefile=paste0('./FBnew.rove_main_age.txt'),
                     intable=table(subdataList$strom_fb@meta.data[,'FB_SubClass02'],subdataList$strom_fb$Location_short), thresh=2)

dir_img="./MP_reScore241104/"
Mono_score=c("FCN1","VCAN","S100A8","S100A12","FYN","ITGAL","CD300E","CCR2","CLU","CX3CR1","FAM65B","LILRB2")  #"LILRB1",CD14","FCGR3A","MARC1","CRISPLD2"
MoDM_score=c("CCL20","C1QA","APOE","PHLDA3","FABP4","CCL18","CD52","BCL22A1A","S100A4","S100A6","S100A11","S100A13","CORO1A","FXDY5","HLA-DRA","HLA-DRB1","S100A8","S100A9","TIMP1","MARCO","TMSB10","HLA-DPB1","HLA-DPA1","HLA-DQA2")
