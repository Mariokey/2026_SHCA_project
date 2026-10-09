# ============================================================================
# Slingshot lineage inference and pseudotime estimation (slingRun)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Trajectory analysis"
# Source     : extracted from Analysis240525.r (L563-754, L4202-4636)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
Intestine_loc=c("ITDM","ITJM","ITIM","ITCM","ITAC","ITTC","ITDC","ITRC")
Stomach_loc=c("GC","GB","GA","GP")
DS_loc=c("Liver","IHB","EHB","GBD","Pancreas","ES")

xxx=table(mtable$System,mtable$Organ)["DS",]
yyy=table(mtable$Organ,mtable$Location_short)[names(xxx[xxx>0]),]
zzz=yyy[,colSums(yyy)>0][c("Liver","IHB","EHB","Gallbladder","Pancreas","Esophagus","Stomach","Intestine"),
                     c(DS_loc, Stomach_loc, Intestine_loc)]

pheatmap(rowPercents(zzz)[,1:dim(zzz)[2]], name = "pct.",
         cluster_rows = F, cluster_cols = F, display_numbers = rowPercents(zzz)[,1:dim(zzz)[2]])

seuloc_Gs=subset(Mobj, Location_short %in% Stomach_loc)
seuloc_ITs=subset(Mobj, Location_short %in% Intestine_loc)

### Stomach
seuloc_Gs=fastRefine(seuloc_Gs, dims = 1:50, nnei = 50, mindist = 0.01, spread = 0.8)
seuloc_Gs$Location_short=droplevels(seuloc_Gs$Location_short)
DimPlot(seuloc_Gs, group.by = "subtypes03", reduction = "umap", label = T, repel = T)
DimPlot(seuloc_Gs, group.by = "seurat_clusters", reduction = "umap", label = T, repel = T)
FeaturePlot(seuloc_Gs, features = "KRT5", reduction = "umap", label = T)
FeaturePlot(seuloc_Gs, features = "PTPRC", reduction = "umap", label = T)

xxx=table(seuloc_Gs$subtypes03)
yyy=100*xxx[xxx>0]/dim(seuloc_Gs)[2]
yyy=sort(yyy,decreasing = T)
tttdf=data.frame(celltype=names(yyy), prop=as.vector(yyy))
tttdf$celltype=factor(tttdf$celltype, levels = unique(tttdf$celltype))
ggplot(tttdf)+geom_bar(aes(celltype, prop, color=celltype, fill=celltype), stat="identity")+
  geom_hline(yintercept = 0.5)+ylab("prop (%)")+theme_cowplot()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
yyy=yyy[yyy<0.5]
DimPlot(subset(seuloc_Gs, subtypes03 %!in% names(yyy)),
        group.by = "subtypes03", reduction = "umap", label = T, repel = T)

seuloc_Gs_Epi=subset(seuloc_Gs, subtypes03 %!in% names(yyy))
seuloc_Gs_Epi$subtypes03=droplevels(seuloc_Gs_Epi$subtypes03)
seuloc_Gs_Epi=subset(seuloc_Gs_Epi, subtypes03 %!in% c("B","DNT","CD4T","CD8T","gdT","NK","NKT","Mono","Mph","Neutro",
                                                         "EC","FB","SMC","MSC","Xendo_unk.","Ximmu_unk."))  #Xepi_unk.
###-------------------------------------------------
seuloc_Gs_Epi=fastRefine(seuloc_Gs_Epi, dims = 1:50, nnei = 80, spread = 1.2)
xxx=DimPlot(seuloc_Gs_Epi, group.by = "subtypes03",
            reduction = "umap", label = T, repel = T, raster = T, raster.dpi = c(300,300))
yyy=DimPlot(seuloc_Gs_Epi, group.by = "seurat_clusters",
            reduction = "umap", label = T, repel = T, raster = T, raster.dpi = c(300,300))
seuloc_Gs_Epi$subtypes03=droplevels(seuloc_Gs_Epi$subtypes03)

fastpHeat(
  table(seuloc_Gs_Epi$seurat_clusters,seuloc_Gs_Epi$Location_short),rowcol="c")

fastpHeat(
  table(seuloc_Gs_Epi$subtypes03,seuloc_Gs_Epi$Location_short),rowcol="c")

FeaturePlot(seuloc_Gs_Epi, features = c("LGR5","CD44","AQP5","LEFTY1", "SMOC2", "ASCL2", "TNFRSF19", "EPHB2", "PROM1","CLDN4"),
            reduction = "umap", label = T, ncol = 2, order = T)

### Intestine
seuloc_ITs=fastRefine(seuloc_ITs, dims = 1:50, nnei = 80, spread = 1.2)
seuloc_ITs$Location_short=droplevels(seuloc_ITs$Location_short)
DimPlot(seuloc_ITs, group.by = "subtypes03", reduction = "umap", label = T, repel = T)
DimPlot(seuloc_ITs, group.by = "seurat_clusters", reduction = "umap", label = T, repel = T)
FeaturePlot(seuloc_ITs, features = "KRT5", reduction = "umap", label = T)
FeaturePlot(seuloc_ITs, features = "PTPRC", reduction = "umap", label = T)
FeaturePlot(seuloc_ITs, features = "PECAM1", reduction = "umap", label = T)
FeaturePlot(seuloc_ITs, features = "LUM", reduction = "umap", label = T)
FeaturePlot(seuloc_ITs, features = c("doublet_scores","pct_counts_mt","nFeature_RNA"), reduction = "umap", label = T)

seuloc_ITs_Epi=subset(seuloc_ITs,
                      subtypes03 %!in% c("B","DNT","CD4T","CD8T","gdT","ILC3",
                                         "NK","NKT","Mono","Mph","Neutro",
                                         "EC","FB","MSC","SMC","Pericytes","Xendo_unk.",
                                         "Xepi_unk.","Ximmu_unk.","Xstrom_unk."))
seuloc_ITs_Epi$subtypes03=droplevels(seuloc_ITs_Epi$subtypes03)
#xxx=FindMarkers(seuloc_ITs_Epi, ident.1 = c("Xendo_unk."), group.by = "subtypes03", test.use = "t", only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)
#xxx=FindMarkers(seuloc_ITs_Epi, ident.1 = c("Xepi_unk."), group.by = "subtypes03", test.use = "t", only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)
#xxx=FindMarkers(seuloc_ITs_Epi, ident.1 = c("unclass.Epi"), group.by = "subtypes03", test.use = "t", only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)
###-------------------------------------------------
seuloc_ITs_Epi=fastRefine(seuloc_ITs_Epi, dims = 1:50, nnei = 80, spread = 1.2)
xxx=DimPlot(seuloc_ITs_Epi, group.by = "subtypes03",
            reduction = "umap", label = T, repel = T, raster = T, raster.dpi = c(300,300))
yyy=DimPlot(seuloc_ITs_Epi, group.by = "seurat_clusters",
            reduction = "umap", label = T, repel = T, raster = T, raster.dpi = c(300,300))
seuloc_ITs_Epi$subtypes03=droplevels(seuloc_ITs_Epi$subtypes03)

fastpHeat(
  table(seuloc_ITs_Epi$seurat_clusters,seuloc_ITs_Epi$Location_short),rowcol="c")

fastpHeat(
  table(seuloc_ITs_Epi$subtypes03,seuloc_ITs_Epi$Location_short),rowcol="c")

FeaturePlot(seuloc_ITs_Epi, features = c("LGR5","ASCL2","SMOC2","RGMB","MKI67","TOP2A","UBE2C","STMN1"), reduction = "umap", label = T, ncol = 2, order = T)

seuloc_ITs_Epi_trajec=subset(seuloc_ITs_Epi, subtypes03 %!in% c("Goblet.Epi","Tuft.Epi"))
seuloc_ITs_Epi_trajec$seurat_clusters=droplevels(seuloc_ITs_Epi_trajec$seurat_clusters)
seuloc_ITs_Epi_trajec@assays$RNA$counts=round(10*seuloc_ITs_Epi_trajec@assays$RNA$counts)
mn2.seuloc_ITs_Epi=runMonocle2(seuloc_ITs_Epi_trajec, "seurat_clusters", mygroup = "Location_short",
                        rmvIdent=c("0","9"), savefix="cds_ITs_Epi_0606.")
tmp_mn2=readRDS(paste0(OutputDir,"MP0605.","monocle_cds_list.dispdeg.rds"))
oplist=list()
oplist=plotMonocle2(tmp_mn2$deg, "Mtype02", group="Sample")
iplist=list()
iplist=plotMonocle2(tmp_mn2$disp, "Mtype02", group="Sample")
MP_mn2_plot01=list(disp=iplist, deg=oplist)

seuloc_Gs_Epi_trajec=subset(seuloc_Gs_Epi, subtypes03 %!in% c("Suprabasal.Epi","Xepi_unk.","Enterocytes"))
seuloc_Gs_Epi_trajec$seurat_clusters=droplevels(seuloc_Gs_Epi_trajec$seurat_clusters)
seuloc_Gs_Epi_trajec@assays$RNA$counts=round(10*seuloc_Gs_Epi_trajec@assays$RNA$counts)
mn2.seuloc_Gs_Epi=runMonocle2(seuloc_Gs_Epi_trajec, "seurat_clusters", mygroup = "Location_short",
                               rmvIdent=c("6","7","8","12"), savefix="cds_Gs_Epi_0606.")

#xxx=FindMarkers(seuloc_Gs_Epi, ident.1 = c("Xepi_unk."), group.by = "subtypes03", test.use = "t", only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)

seuloc_ITs_Epi_trajec=subset(seuloc_ITs_Epi_trajec, seurat_clusters %!in% c(0,9))
seuloc_ITs_Epi_trajec=RunUMAP(seuloc_ITs_Epi_trajec,dims = 1:50, reduction="pca_harmony", reduction.name = "umap",
                             n.neighbors = 80, min.dist = 0.15, spread = 1.2, verbose = FALSE)
DimPlot(seuloc_ITs_Epi_trajec,
  group.by = "subtypes03", label = T, raster = F, repel = T)+NoTitle

seuloc_Gs_Epi_trajec=subset(seuloc_Gs_Epi_trajec, seurat_clusters %!in% c(6,7,8,12))
seuloc_Gs_Epi_trajec=RunUMAP(seuloc_Gs_Epi_trajec,dims = 1:50, reduction="pca_harmony", reduction.name = "umap",
                             n.neighbors = 80, min.dist = 0.15, spread = 1.2, verbose = FALSE)
DimPlot(seuloc_Gs_Epi_trajec,
  group.by = "subtypes03", label = T, raster = F, repel = T)+NoTitle

plotobj=mn2.seuloc_Gs_Epi$disp

state_order=plotobj$heat_comb@matrix %>% rownames()
color_state=hue_pal()(length(state_order))
names(color_state)=state_order
ct_order=plotobj$heat_comb@matrix %>% colnames()
color_tmp=hue_pal()(length(ct_order))
names(color_tmp)=ct_order

#plotobj$heat_comb
#plotobj$state
tmpxxx=plotobj$pseudotime
tmpxxx$data$Pseudotime=max(tmpxxx$data$Pseudotime)-tmpxxx$data$Pseudotime

tmpxxx=plotobj$split_by_phase
tmpxxx$facet$params$nrow=2
tmpxxx+myfont+NoAxes()

tmpxxx=plotobj$counts_split_by_phase_celltype_and_state
tmpDF=data.frame()
group=c("GC", "GB", "GA", "GP")
tmpDF=lapply(group,FUN = function(x){tmpxxx[,,x] %>% rowSums()}) %>%
  as.data.frame.list() %>% set_colnames(.,group)
ggplot(colPercents(tmpDF)[1:dim(tmpDF)[1],] %>% melt() %>% set_colnames(.,c("State","Group","Fraction")))+
  geom_bar(aes(Group, Fraction, fill=as.factor(State), color=as.factor(State)), stat = "identity", position = position_stack())+
  scale_fill_manual(values = color_state)+scale_color_manual(values = color_state)+
  xlab("")+theme_minimal()+myfont+NoLegend()

tmpxxx=plotobj$celltype$data
tmpxxx$celltypes=tmpxxx$seurat_clusters
time_mid=tmpxxx %>% group_by(celltypes) %>% summarise_at("Pseudotime",mean)
tmpxxx$celltypes=factor(tmpxxx$celltypes,levels = time_mid[order(time_mid$Pseudotime),]$celltypes)
ggplot(tmpxxx)+
  geom_violin(aes(celltypes, Pseudotime),fill="lightgrey", color="lightgrey", scale = "width", trim = FALSE, adjust=2)+
  geom_boxplot(aes(celltypes, Pseudotime, fill=celltypes), width=0.3, outlier.alpha = 0, outlier.size = 0)+
  scale_fill_manual(values = color_tmp)+theme_minimal()+NoLegend()+xlab("")+coord_flip()+myfont

tmpxxx=plotobj$counts_split_by_celltype_and_state
pheatmap(rowPercents(tmpxxx)[,1:dim(tmpxxx)[2]], cluster_rows = F, cluster_cols = F,
         display_numbers = rowPercents(tmpxxx)[,1:dim(tmpxxx)[2]], name="Frac.")

### slingshot analysis 0704
Idents(seuloc_Gs_Epi_trajec)="seurat_clusters"
seuloc_Gs_Epi_trajec$Cellname=rownames(seuloc_Gs_Epi_trajec@meta.data)
seuloc_Gs_Epi_trajec$Location_short=droplevels(seuloc_Gs_Epi_trajec$Location_short)
GsEpi_sls_res0704=slingRun(seuloc_Gs_Epi_trajec, samplesize=1, CTname="seurat_clusters", cellid="Cellname",
                        start.type=c("9"),
                        mycolor89=NULL, savesce=paste0(resDir,"/GsEpi.slingshot.0704.rds"))
sce=readRDS(paste0(resDir,"/GsEpi.slingshot.0704.rds"))
GsEpi_tradeResList=tradeS(sce, samplesize=1, SvE=FALSE,  Asso=TRUE, outdirprefix=paste0(resDir,"/tradeSeq0704/GsEpi"))
save.image("./data/R_data/PRJ21_20240525.Rdata")

Idents(seuloc_ITs_Epi_trajec)="seurat_clusters"
seuloc_ITs_Epi_trajec$Cellname=rownames(seuloc_ITs_Epi_trajec@meta.data)
seuloc_ITs_Epi_trajec$Location_short=droplevels(seuloc_ITs_Epi_trajec$Location_short)
ITsEpi_sls_res0704=slingRun(seuloc_ITs_Epi_trajec, samplesize=1, CTname="seurat_clusters", cellid="Cellname",
                           start.type=c("10"),
                           mycolor89=NULL, savesce=paste0(resDir,"/ITsEpi.slingshot.0704.rds"))
sce=readRDS(paste0(resDir,"/ITsEpi.slingshot.0704.rds"))
ITsEpi_tradeResList=tradeS(sce, samplesize=1, SvE=FALSE,  Asso=TRUE, outdirprefix=paste0(resDir,"/tradeSeq0704/ITsEpi"))

name="regulation of transmembrane transporter activity" #"chromosome segregation" #"cell cycle process" #"lipid transport"
name_go="GO:0022898"
JAKplot=gseaplot(GsEpi_tradeResList$L1_GSEres$gse_obj, by = "all", title = name, geneSetID = name_go)
addDFtop10=head(JAKplot$patches$plots[[1]]$data[JAKplot$patches$plots[[1]]$data$ymax!=0,], 10)
genes_show=names(GsEpi_tradeResList$L1_GSEres$gse_obj@geneList)[addDFtop10$x]
addDFtop10$gene=names(genes_show)
JAKplot$patches$plots[[1]]=JAKplot$patches$plots[[1]]+
  ggrepel::geom_text_repel(data=addDFtop10, aes(x=x, y=ymin, label=gene), direction="x", color="red",
                           force = 1, force_pull = 0, max.overlaps = 1000, min.segment.length = 0.1, nudge_x = 0, nudge_y = -0.2, angle=90, hjust=1)

figfix="_ITs"
tmpTrade=ITsEpi_tradeResList

GSEdir=paste0(resDir,"/tradeSeq0704/","GSE_pathway.Early5Late5",figfix,"/")
dir.create(GSEdir, recursive = T)
LLLnames=names(tmpTrade)[grepl("_smoothExp$",names(tmpTrade))]
LLLnames=gsub("_smoothExp$","",LLLnames)
for(Lidx in 1:length(LLLnames)){
  #  next
  #}
  LgseRes=paste0(LLLnames[Lidx],"_GSEres")
  LgseSmooth=paste0(LLLnames[Lidx],"_smoothExp")
  original_gene_list=rownames(tmpTrade[[LgseSmooth]])
  gseDF=tmpTrade[[LgseRes]]$gse_fdf
  gseOBJ=tmpTrade[[LgseRes]]$gse_obj
write.xlsx(tmpTrade[[LgseSmooth]][,1:3] %>% as.data.frame(), rowNames=TRUE, file = paste0(GSEdir,"/L",Lidx,".smoothHeatmap_col123.xlsx"))
}

### slinshot functions

gseFun=function(gene_list,...){

  ## omit any NA values
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list,
               ont ="BP",
               keyType = "SYMBOL",
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = "org.Hs.eg.db",
               pAdjustMethod = "none")

  gseDf=gse@result
  gseDf=gseDf[gseDf$pvalue<=0.05,]

  return(
    list(
      gse_obj=gse,
      gse_fdf=gseDf
    )
  )
}

slingRun=function(tmpData, samplesize=1, CTname="CT89", cellid="Cellname", start.type=c("Basal_IGFBP3"), end.type=NULL, mycolor89=NULL,
                  savesce="DATA_DIR", ...){

  #### slingshot 20230923
  library(slingshot)
  library(SingleCellExperiment)

  tmpData=subset(tmpData, cells = sample(colnames(tmpData), size = round(samplesize*ncol(tmpData))))
  sls_res_list=list()

  Idents(tmpData)=CTname
  sce <- as.SingleCellExperiment(tmpData)

  sce <- slingshot(sce, clusterLabels = "ident", reducedDim = "UMAP", approx_points = 200, reweight = FALSE, reassign = FALSE,
                   dist.method = "slingshot", maxit=20, stretch=0, thresh=50, shrink=1, shrink.method="cosine", omega=TRUE, omega_scale=3,
                   allow.breaks = FALSE, extend = 'n', start.clus=start.type, end.clus = end.type)
  if(samplesize!=1){
    savesce_seu=gsub(".rds$",".subseu.rds",savesce)
    saveRDS(tmpData, savesce_seu)
  }
  saveRDS(sce, savesce)

  print("slingshot done, start plotting...")
  slingshot_df <- colData(sce)
  col_n=colnames(slingshot_df)
  col_n=col_n[grepl("^slingP",col_n)]
  c_slingshot_df=data.frame(Cellname=rownames(slingshot_df))   # cellid NAME must be "Cellname"
  for(x in col_n){
    slingshot_df_tmp = data.frame(row.names = rownames(slingshot_df),
                                  value=as.vector(unlist(slingshot_df[x])))
    colnames(slingshot_df_tmp)=x
    c_slingshot_df=cbind(c_slingshot_df, slingshot_df_tmp)
  }

  plotdf=as.data.frame(reducedDims(sce)$UMAP)
  plotdf$Cellname=rownames(plotdf)
  metadf=tmpData@meta.data[, c(cellid,"Organ",CTname)]
  print(head(plotdf))
  print(head(c_slingshot_df))
  plotdf=left_join(plotdf, metadf, by=cellid)
  plotdf=left_join(plotdf, c_slingshot_df, by=cellid)

  sls_res_list[["plot_data"]]=plotdf

  if(is.null(mycolor89)){
    mycolor89=sample(scPioneer::scPalette2(length(levels(tmpData@meta.data[, CTname]))))
    names(mycolor89)=levels(tmpData@meta.data[, CTname])
  }

  mycurves=slingCurves(sce)
  curvAll=data.frame()
  xxx=names(mycurves)
  for(i in 1:length(xxx)){
    tmpcurv=as.data.frame(mycurves[xxx[i]][[1]]$s)
    tmpcurv$Lineage=rep(i, dim(tmpcurv)[1])
    if(i==1){
      curvAll=tmpcurv
    }else{
      curvAll=rbind(curvAll,tmpcurv)
    }
  }
  sls_res_list[["plot_curve_data"]]=curvAll

  AAA=ggplot(plotdf)+geom_point(aes_string("umap_1", "umap_2", color=CTname), size=0.2)+
    scale_color_manual(values = mycolor89)+#guides(color = guide_legend(title = "Celltype", override.aes = list(size=5)))+
    geom_path(data = curvAll, aes(umap_1, umap_2, size=as.factor(Lineage)), color="black", arrow = arrow(type="closed", angle = 15, length = unit(0.1, "inches")))+
    scale_size_manual(name="Lineage", breaks = sort(unique(curvAll$Lineage)), values = c(1.5,1,0.75,0.5))+
    guides(color = guide_legend(title = "Celltype", override.aes = list(size=5), order = 1),
           size = guide_legend(title = "Lineage", override.aes = list(), order = 2))+
    theme_minimal_grid()
  sls_res_list[["UMAP_allCurve"]]=AAA

  for(i in sort(unique(curvAll$Lineage))){
    plotdf$Pseudotime=as.numeric(rescale(plotdf[,paste0("slingPseudotime_",i)], to = c(0,100)))
    plotdf$Lineage=rep(paste0("Lineage", i), dim(plotdf)[1])
    p2=ggplot(plotdf)+geom_density(aes(Pseudotime, color=Organ), adjust=2)+theme_cowplot()+theme(legend.position = "bottom")
    zzz=ggplot(plotdf)+geom_point(aes_string("Lineage","Pseudotime", color=CTname), position = position_jitter())+
      scale_color_manual(values = mycolor89)+guides(color = guide_legend(title = "Celltype"))+
      theme_minimal_grid()+coord_flip()+theme(legend.position = "top", axis.title = element_blank(), axis.text = element_blank())+NoLegend()
    mmm0=plot_grid(zzz,p2, nrow = 2,rel_heights = c(0.2,1), align = "v")
    sls_res_list[[paste0("densityPoint_Curve", i)]]=mmm0
  }

  return(sls_res_list)

}

tradeS=function(sce, samplesize=0.05, SvE=FALSE,  Asso=TRUE, outdirprefix="DATA_DIR", ...){

  ### Tradeseq analysis : within 30min < 4000 cells
  library(tradeSeq)
  library(RColorBrewer)
  library(SingleCellExperiment)
  library(slingshot)
  library(scran)
  library(org.Hs.eg.db)

  tradeResList=list()

BPPARAM = BiocParallel::bpparam()
BPPARAM$workers=6

sce1k=sce[, sample(colnames(sce), size=round(samplesize*length(colnames(sce))))]

tseq_counts=sce1k@assays@data@listData$counts
tseq_pseudotime <- slingPseudotime(sce, na = FALSE)  ####!!! na=FALSE: arclength along each curve will be returned for NA cells.
if(dim(tseq_pseudotime)[2]==1){
  tseq_pseudotime=data.frame(Lineage1=as.vector(tseq_pseudotime[colnames(sce1k),]),
                             row.names = names(tseq_pseudotime[colnames(sce1k),]))
}else{
  tseq_pseudotime=tseq_pseudotime[colnames(sce1k),]
}
tseq_cellWeights <- slingCurveWeights(sce1k)

tradeResList[["tradeS_sce1k"]]=sce1k
tradeResList[["tradeS_sce1k_time"]]=tseq_pseudotime

set.seed(7)
x <- org.Hs.egGENETYPE
mapped_type <- mappedkeys(x)
xx <- as.list(x[mapped_type])
pcgene_id=names(xx[xx=="protein-coding"])
y <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(y)
yy <- as.list(y[mapped_genes])
pcgene_symbol=as.vector(unlist(yy[pcgene_id]))

variable_genes <-
  sce1k %>%
  scran::modelGeneVar() %>%
  scran::getTopHVGs(prop=0.1)

tseq_selectGenes=variable_genes  ### or from obj seurat data
tseq_selectGenes=tseq_selectGenes[tseq_selectGenes %in% pcgene_symbol]
tseq_sce <- fitGAM(counts = tseq_counts, pseudotime = tseq_pseudotime, cellWeights = tseq_cellWeights, #conditions = factor(obj$Group),
                   parallel = TRUE, BPPARAM = BPPARAM, control = mgcv::gam.control(),
                   genes=tseq_selectGenes, nknots = 7, verbose = T)
tradeResList[["tradeS_sce1k_tseq_sce"]]=tseq_sce

num_Vgene_used=length(tseq_selectGenes)
print(paste0("Number of input variable genes for tradeS: ", num_Vgene_used))
table(rowData(tseq_sce)$tradeSeq$converged)

if(Asso){
  assoRes <- associationTest(tseq_sce, lineages=TRUE)
  tradeResList[["tradeS_assoRes"]]=assoRes
}

if(SvE){
  startRes <- startVsEndTest(tseq_sce, lineages=TRUE)  ##pseudotimeValues = c(0, 1)  ##for large data l2fc = log2(2) to reduce FalsePositiveRate
  tradeResList[["tradeS_startRes"]]=startRes
}

  if(Asso){
    tempAll=as.data.frame(t(tseq_counts[rownames(assoRes),rownames(tseq_pseudotime)]))
    tempAll$cellid=rownames(tempAll)

    for(mark in 1:dim(tseq_pseudotime)[2]){

      targetL=assoRes[assoRes[,paste0("pvalue_",mark)]<0.05,]
      temp=data.frame(cellid=rownames(tseq_pseudotime), #names(selectTime),
                      Lineage=as.vector(tseq_pseudotime[,paste0("Lineage",mark)]))  #as.vector(selectTime))
      intersctGenes=intersect(colnames(tempAll),rownames(targetL))
      temp=temp[order(temp[,"Lineage"], decreasing = F),]
      temp=left_join(temp, tempAll[,c("cellid", intersctGenes)], by="cellid")

      tradeResList[[paste0("L",mark,"_rawExp")]]=temp

      nana=data.frame(Lineage=temp$Lineage)
      for(i in 3:dim(temp)[2]){
        temp$gene=temp[,i]
        loessmod <- loess(gene ~ Lineage, data = temp, span = 1)
        smoothed_data <- predict(loessmod)
        smooth_scale <- rescale(smoothed_data, to=c(0,100))
        nana <- cbind(nana, smooth_scale)
      }
      temp$gene=NULL
      temp$cellid=NULL
      colnames(nana)=colnames(temp)
      #tradeResList[[paste0("L",mark,"_smoothExp")]]=nana
tmpheat=t(as.matrix(nana[,-1]))
colnames(tmpheat)=nana$Lineage
gene_weights=c()
pseudotime=as.numeric(colnames(tmpheat))
for(i in 1:dim(tmpheat)[1]){
  gene_val=(tmpheat[i,]+1)/sum(tmpheat[i,])
  gene_weights=c(gene_weights, sum(gene_val*pseudotime))
}
names(gene_weights)=rownames(tmpheat)
gene_weights=sort(gene_weights,decreasing = T)
AAA=pheatmap(tmpheat[names(gene_weights),], color = viridisLite::viridis(256, option = "B"), name = "Relative\nexpression",
             heatmap_legend_param = list(direction = "horizontal", title_position = "leftcenter"),
             cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
AAA=draw(AAA, heatmap_legend_side = "bottom")

tradeResList[[paste0("L",mark,"_smoothExp")]]=tmpheat[names(gene_weights),]
tradeResList[[paste0("L",mark,"_smoothExp_heatplot")]]=AAA

      original_gene_list=gene_weights   ###avg_log2FC MUST decreasing=TRUE
      original_gene_list=rescale(original_gene_list, to=c(-1,1))
      res.GSE=gseFun(original_gene_list)
      res.GSE$gse_fdf=res.GSE$gse_fdf[order(res.GSE$gse_fdf$NES, decreasing = T),]  ### !!!! IMPORTANT, NES>0 are enrichment of genes at late pseudotime
      head(res.GSE$gse_fdf)
      tradeResList[[paste0("L",mark,"_GSEres")]]=res.GSE
    }
  }

  saveRDS(tradeResList, file = paste0(outdirprefix,".tradeseqResult.rds"))
  return(tradeResList)
}
