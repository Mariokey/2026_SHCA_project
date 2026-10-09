# ============================================================================
# Monocle2 trajectory inference (runMonocle2)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Trajectory analysis"
# Source     : extracted from Analysis240525.r (L1663-1819, L5883-5993)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
runMonocle2=function(obj, myident, mygroup="group", rmvIdent=NULL, savefix="tmp.", mode=c("deg","disp"), ...){

library(monocle)

DefaultAssay(obj)="RNA"
obj=obj[,obj@meta.data[,myident] %!in% rmvIdent]
obj@meta.data[,myident]=droplevels(obj@meta.data[,myident])

#Extract data, phenotype data, and feature data from the SeuratObject
monodata <- as(as.matrix(obj@assays$RNA$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = obj@meta.data)
fData <- data.frame(gene_short_name = row.names(monodata), row.names = row.names(monodata))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(monodata,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 10)
expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))

saveRDS(object = monocle_cds, file = paste0(savefix,"monocle_cds.rds"))

### Trajectory step 1: HVG or DEG (optional?), or provided!
if("disp" %in% mode){
  monocle_cds_disp=monocle_cds
  disp_table <- dispersionTable(monocle_cds_disp)
  unsup_clustering_genes1 <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1.5 * dispersion_fit)
  monocle_cds_disp <- setOrderingFilter(monocle_cds_disp, unsup_clustering_genes1$gene_id)
  monocle_cds_disp <- reduceDimension(monocle_cds_disp, max_components = 2,
                                      method = 'DDRTree')
  monocle_cds_disp <- orderCells(monocle_cds_disp)
}else{monocle_cds_disp=NULL}
if("deg" %in% mode){
  monocle_cds_deg=monocle_cds
  diff_test_res <- differentialGeneTest(monocle_cds_deg[expressed_genes,],
                                        fullModelFormulaStr = paste0("~",myident))   ###???percent.mt
  diff_ordering_genes <- row.names (subset(diff_test_res, qval < 0.01 & num_cells_expressed > 30))
  monocle_cds_deg <- setOrderingFilter(monocle_cds_deg, diff_ordering_genes)
  monocle_cds_deg <- reduceDimension(monocle_cds_deg, max_components = 2,
                                     method = 'DDRTree')
  monocle_cds_deg <- orderCells(monocle_cds_deg)
}else{monocle_cds_deg=NULL}
calc_cds=list(disp=monocle_cds_disp, deg=monocle_cds_deg)
saveRDS(object = calc_cds, file = paste0(savefix,"monocle_cds_list.dispdeg.rds"))

plotMonocle2=function(xxx, myident, group="group", ...){

  plist=list()

  #xxx=monocle_cds_deg
  pData(xxx)[,myident]=droplevels(pData(xxx)[,myident])
  pData(xxx)[,group]=droplevels(pData(xxx)[,group])

  mylen=length(levels(pData(xxx)[,group]))
  plist[["split_by_phase"]]=plot_cell_trajectory(xxx, color_by = myident)+theme(legend.position = "none")+facet_wrap(as.formula(paste("~", group)), nrow = mylen)    ### +facet_wrap(enquote(group_var), nrow = 3)
  plist[["split_by_phase_and_celltype"]]=plot_cell_trajectory(xxx, color_by = myident)+theme(legend.position = "none")+facet_wrap(as.formula(paste("~", group, "+", myident)), nrow = mylen)

  plist[["celltype"]]=plot_cell_trajectory(xxx, color_by = myident)
  plist[["split_by_celltype"]]=plot_cell_trajectory(xxx, color_by = myident)+
    facet_wrap(as.formula(paste("~", myident)), nrow = if_else(length(unique(pData(xxx)[,myident]))<=3, 1, floor(length(unique(pData(xxx)[,myident]))/3)+1))+NoLegend()
  plist[["state"]]=plot_cell_trajectory(xxx, color_by = "State")
  plist[["pseudotime"]]=plot_cell_trajectory(xxx, color_by = "Pseudotime")

  tttdf=table(pData(xxx)$State, pData(xxx)[,myident])
  plist[["counts_split_by_celltype_and_state"]]=tttdf
  plist[["counts_split_by_phase_celltype_and_state"]]=table(pData(xxx)$State, pData(xxx)[,myident], pData(xxx)[,group])

  plist[["heat_comb"]]=pheatmap(colPercents(plist[[7]])[1:dim(plist[[7]])[1],], cluster_rows = F, cluster_cols = F,
                                display_numbers = colPercents(plist[[7]])[1:dim(plist[[7]])[1],], name="Frac.")

  mdf=as.data.frame(plist[[8]])
  colnames(mdf)=c("State","Celltype","Group","Counts")
  myphase=levels(pData(xxx)[,group]) #unique(mdf$Group)
  for (j in 1:length(myphase)){
    #j=1
    pa=myphase[j]
    tdf=mdf[mdf$Group==pa,]

    tmat=as.data.frame(tidyr::pivot_wider(tdf, id_cols="State", names_from = "Celltype", values_from = "Counts"))
    rownames(tmat)=tmat$State
    tmat=tmat[,-1]

    #heatlist[[paste0(pa,"_counts")]]=gridExtra::grid.table(tmat)
    plist[[paste0("heat_",pa)]]=pheatmap(colPercents(tmat)[1:dim(tmat)[1],], cluster_rows = F, cluster_cols = F,
                                         display_numbers = colPercents(tmat)[1:dim(tmat)[1],], name="Frac.")
  }

  return(plist)
}

oplist=list()
if("deg" %in% mode){
oplist=plotMonocle2(monocle_cds_deg, myident, group=mygroup)
}
iplist=list()
if("disp" %in% mode){
iplist=plotMonocle2(monocle_cds_disp, myident, group=mygroup)
}
plist=list(disp=iplist, deg=oplist)
saveRDS(object = plist, file = paste0(savefix,"monocle_cds_plotlist.dispdeg.rds"))

return(plist)

}

### trajectory 1023
dir.create("./MP_mn2/")
mp_trajec=subdataList$myl_mp
mp_trajec=subset(mp_trajec, MYL_SubClass04_bak %!in% c("rmv_BS","rmv_Intestine","rmv_skin","rmv_spleen","rmv_UDT"))
DimPlot(mp_trajec, group.by = "MYL_SubClass04_bak", label = T, repel = T,raster = T)+NoLegend()
mp_trajec=RunUMAP(mp_trajec,dims = 1:20, reduction="pca_harmony", reduction.name = "umap",
                  n.neighbors = 150, min.dist = 0.35, spread = 1.2, verbose = FALSE)
data_list=list()
data_list[["DS1"]]=c("Liver","BS")
data_list[["DS2"]]=c("UDT","Intestine")
data_list[["RS"]]=c("Lung","Airway")
data_list[["US"]]=c("UB","Kidney")
data_list[["CS"]]=c("Spleen","Artery")
for(org_idx in names(data_list)){
  xxx_organ=data_list[[org_idx]]
  mp_sub_trajec=subset(mp_trajec, Organ %in% xxx_organ)
  mp_sub_trajec=subset(mp_sub_trajec, cells = sample(colnames(mp_sub_trajec), size = round(0.3*ncol(mp_sub_trajec))))
mp_sub_trajec$MYL_SubClass04_bak=droplevels(mp_sub_trajec$MYL_SubClass04_bak)
mp_sub_trajec$Organ=droplevels(mp_sub_trajec$Organ)
mp_sub_trajec@assays$RNA$counts=round(10*mp_sub_trajec@assays$RNA$counts)
print(xxx_organ)
print(mp_sub_trajec)
mn2.mp_sub_trajec=runMonocle2(mp_sub_trajec, "MYL_SubClass04_bak", mygroup = "Organ", mode = c("disp"),
                              rmvIdent=NULL, savefix=paste0("./MP_mn2/cds_mp_trajec_1023.",org_idx,"."))
mn2_organ_list[[paste(xxx_organ,collapse = "",sep="")]]=mn2.mp_sub_trajec
}

mp_trajec=reDRClusterTest(mp_trajec, sp="hum", ndim=11, nClu=120, mdist=0.15,
                          formalFeature=NULL, mypca="harmony", harmony=0, prepro=FALSE,
                          mybatch="sample", umap=TRUE, vstN=2000, regVar=c("pct_counts_mt"), rmvHVGribo=TRUE)
mp_trajec=reResClu(mp_trajec, reClu = TRUE, res = c(0.8,1.2), snnres = "RNA_snn_res.0.8")
mp_trajec$MYL_SubClass05=mp_trajec$MYL_SubClass04_bak
mp_trajec$MYL_SubClass05=droplevels(mp_trajec$MYL_SubClass05)
mp_trajec$Organ=droplevels(mp_trajec$Organ)
oldID=levels(mp_trajec$MYL_SubClass05)
newID=c("Mph04_ITGAX","DC01_CD1C","Mph02_C1QB","Mph01_STARD13","Mph02_C1QB",
        "DC02_CD207","Mono01_VCAN","Mono01_VCAN","Mono01_VCAN","Mono01_VCAN",
        "Mph02_C1QB","Mph01_STARD13","Mph02_C1QB","Mph02_C1QB","Mph02_C1QB",
        "Mph02_C1QB","Mph02_C1QB","Mph02_C1QB","Mph02_C1QB","Mph02_C1QB",
        "Mph02_C1QB","Mph02_C1QB","Mph02_C1QB","Mph02_C1QB","Mph01_STARD13")
for(x in 1:length(oldID)){
mp_trajec=modifyAno(mp_trajec, anoCol = "MYL_SubClass05", oldID = oldID[x], newID = newID[x] )
}
DimPlot(mp_trajec, group.by = "MYL_SubClass05", label = T, raster = F)
mp_trajec=reResClu(mp_trajec, reClu = TRUE, res = c(1), snnres = "RNA_snn_res.1")
DimPlot(mp_trajec, group.by = "seurat_clusters", label = T, raster = F)
DimPlot_idx(mp_trajec, group.by = "seurat_clusters", cols=scPalette2(30))
DimPlot_idx(mp_trajec, group.by = "MYL_SubClass04_bak", cols=scPalette2(30))

DEG2FUN_mp_trajec_clusters=DEG2FUN(mp_trajec, myiden="seurat_clusters", groupname="Location_short",
                              pval_threshold = 0.01, logfc_threshold = 0.25, min_pct = 0.1, DEGenefilter=FALSE,
                              DEGct=TRUE, DEGgroup=FALSE,
                              maxCellperIdents=2000,
                              FunctionGO=TRUE, FunctionKEGG=TRUE, FunctionGSE=TRUE,
                              OutDir = "./MP_mn2/", suffix = paste0("mp_trajec","_subcts"))

### trajectory 1104
dir.create("./MP_reScore241104/")
FinalFig="./MP_reScore241104/"
MonoDP_gs=c(#"PRSS57","CDK6", "HSPD1", "LGALS1","CEBPA",
            "FCN1", "VCAN", "CD14", "FCGR3A", "S100A4", "S100A6", "S100A8", #"CD68",
            "CST3", "SCT","TGFBI","SAMHD1", "FCER1A","CPVL"#,
            #"CENPK","MYBL2", "RHAG","SLC14A1", "GINS2","MCM10","AURKB","CDC45"
            )
mp_trajec=AddModuleScore(mp_trajec, features = list(MonoDP_gs),
                     ctrl = 100, name = "MonoDP0")
FeaturePlot2(mp_trajec, features = c("MPresd01","MonoDP01"), reduction = "umap",label = F, ncol = 2)

mp_trajec_cluDegs=FindAllMarkers(mp_trajec, group.by = "seurat_clusters", test.use = "t", logfc.threshold = 0.25,
                   max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
Idents(mp_trajec)="MYL_SubClass04_bak"
mp_trajec_class04Degs=FindAllMarkers(mp_trajec, group.by = "MYL_SubClass04_bak", test.use = "t", logfc.threshold = 0.25,
                                 max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)

mp_trajec=modifyAno(mp_trajec, anoCol = "seurat_clusters", oldID = 17, newID = 0)
mp_trajec$seurat_clusters=factor(mp_trajec$seurat_clusters, levels = 0:15)
mp_trajec=reDRClusterTest(mp_trajec, sp="hum", ndim=13, nClu=150, mdist=0.2,
                          formalFeature=NULL, mypca="harmony", harmony=0, prepro=FALSE,
                          mybatch="sample", umap=TRUE, vstN=2000, regVar=c("pct_counts_mt"), rmvHVGribo=TRUE)

dfin_x=mp_trajec_cluDegs
dfin_x$group=rep("group1",dim(dfin_x)[1])
dfin_y=mp_trajec_class04Degs
dfin_y$group=rep("group2",dim(dfin_y)[1])
dfin2=rbind(dfin_x,dfin_y)
dfin2=dfin2[dfin2$p_val_adj<=0.01
            & dfin2$avg_log2FC>=0.35
            & (dfin2$pct.1-dfin2$pct.2)>=0.15
            ,]
JACdfin=CorTable(dfin2, JC = TRUE)
pheatmap(JACdfin, cluster_rows = F, cluster_cols = F, display_numbers = JACdfin, name = "JacDis.")

Idents(mp_trajec)="MYL_SubClass04_bak"
aveScore1=mp_trajec@meta.data[,c("seurat_clusters","MYL_SubClass04_bak","MPresd01","MonoDP01")]
ggplot(aveScore1)+geom_point(aes(MPresd01, MonoDP01, color=seurat_clusters, fill=seurat_clusters), alpha=0.2, shape=22)+
  facet_wrap(~seurat_clusters, ncol = 4) + cowplot::theme_cowplot()

xxx=plot2dfun(mp_trajec,
              features=c("MPresd01","MonoDP01","MKPP_score1"))
