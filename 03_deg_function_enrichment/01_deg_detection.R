# ============================================================================
# Marker / DEG detection per cluster, per location x subtype, and cluster-unique DEG filtering
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Function annotation"
# Source     : extracted from Analysis240525.r (L508-545, L2521-2600, L5465-5497, L7350-7379, L11120-11145)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
mtable=Mobj@meta.data
mtable=mtable[order(mtable$index_order00,decreasing = FALSE),]
Mobj$Location_short=factor(Mobj$Location_short, levels = as.vector(unique(mtable$Location_short)))
Mobj$System=factor(Mobj$System, levels = as.vector(unique(mtable$System)))

nmtx=table(Mobj$subtypes03,Mobj$Location_short)
rmtx=colPercents(nmtx)[1:dim(nmtx)[1],]
pheatmap(rmtx, cluster_rows = F, cluster_cols = F, display_numbers = rmtx)

Idents(Mobj)="subtypes03"
subtypes03CTdegs_RNA = FindAllMarkers(Mobj, group.by = "subtypes03",
                                      test.use = "t", only.pos = TRUE, assay = "RNA", max.cells.per.ident = 2000)

DEGcl=subtypes03CTdegs_RNA
DEGcl$pct.diff=DEGcl$pct.1-DEGcl$pct.2
DEGcl=DEGcl[!grepl("^ENSG|^LINC|^A[PLC]\\d\\d\\d",DEGcl$gene) & DEGcl$pct.1>0.3,]
DEGcl_top=as.data.frame(DEGcl %>% group_by(cluster) %>% top_n(n=5, wt=(avg_log2FC+5*pct.diff)))
DEGcl_top=DEGcl_top[grepl(".Epi$|Hepatocytes|Enterocytes",DEGcl_top$cluster),]
xxx=DotPlot(Mobj, features = unique(DEGcl_top$gene), group.by = "subtypes03")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
runFig(paste0(OutputDir,"/", "XXX.DotPlot_DEG.epi"), xxx, 28, 10)
DEGcl_top=as.data.frame(DEGcl %>% group_by(cluster) %>% top_n(n=8, wt=(avg_log2FC+5*pct.diff)))
DEGcl_top=DEGcl_top[grepl("_unk.$",DEGcl_top$cluster),]
xxx=DotPlot(Mobj, features = unique(DEGcl_top$gene), group.by = "subtypes03")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
runFig(paste0(OutputDir,"/", "XXX.DotPlot_DEG.X_unk"), xxx, 12, 10)

select_major_markers=c("MS4A1","CD79A","THEMIS","CD3D","CD4","CD8B","TRGC2","TRDC","KLRD1","NKG7","IL4I1",
                       "FCN1","C1QC","CSF3R",
                       "DCT","ADGRL4","PECAM1","NT5E","CLEC3B","COL1A2","DCN","ACTA2","TAGLN","ABCC9","NR2F2",
                       "MAGEB2","TCF3","SPO11","SYCP3",'SPACA1','ACRV1',"TNP1","PRM3","MKI67","TOP2A")
DotPlot(Mobj, features = select_major_markers, group.by = "subtypes03")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")

#xxx=DoHeatmapSelf(subset(Mobj, cells = sample(colnames(Mobj), size = round(0.05*ncol(Mobj)))),

##############################
Idents(Mobj)="refine_ct00"
DEGs_Mobj_33CT=FindAllMarkers(Mobj, group.by = "refine_ct00", only.pos = T, max.cells.per.ident = 2000, test.use = "t")

formal_major_markers=c(
  "MS4A1","MZB1",
  "CD3D","KLRD1",
  "FCN1","C1QC","FCGR3B","CPA3",

  "CDH1","EPCAM","LGR4","KRT18",#"KRT8",
  "KRT19","CD24","ALB","HP", #BEC/Hep  #"CFTR","SOX9"
  "PRSS1","PRSS2",  #Pancreas
  "NEUROD1","RFX6",  #Endocrine
  "KRT7","CLDN4", #Collecting-ductal
  "PGC","PHGR1","PIGR", #Gastric/Intestinal
  "WFDC2", #Tracheal
  "CALML3", #Skin
  "KRT5","KRT15","KRT1","KRT10", #Super/Basal
  "SERPINE2","PRKAR2B", #Granulosa #GSTA1
  "CHAT","POU2F3", #Tuft
  "PRG4","WT1", #Mesothelial

  "EGFL7","COL1A2","ACTA2",
  "TCF3","SYCP3","TNP1",  #TEX15/MAGEB2
  "MKI67"
)
DotPlot2(Mobj,features = formal_major_markers, group.by = "refine_ct00")

### Granulosa->ST
### RPS_Epi_undef4->Granulosa
### Epi_undef2->remove (RPS/VIM/FTL, no epi markers)
### Epi_undef3->Epi_undef2

Mobj=subset(Mobj, refine_ct00 %!in% c("Epi_undef2"))
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Granulosa", newID = "ST")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "RPS_Epi_undef4", newID = "Granulosa")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Epi_undef3", newID = "Epi_undef2")

CT33_order=c("B","Plasma",
             "T","NK",
             "Mono","Mph","Neutro","Mast",
             "BEC","Hepatocyte","Pancreas",
             "Endocrine","Collecting_ductal",
             "Gastric","Intestinal","Tracheal",
             "Skin","Basal","Superbasal",
             "Granulosa","Tuft","Mesothelial",
             "Epi_undef1","Epi_undef2",
             "EC","FB","SMC","SPG","SPC","ST","Prolif." )
Mobj$refine_ct00=factor(Mobj$refine_ct00, levels = CT33_order)

Mobj=subset(Mobj, subtypes01 %!in% c("granulosa","Xgerm_unk.","Xstrom_unk."))
Mobj$subtypes011=Mobj$subtypes01
Mobj=modifyAno(Mobj, anoCol = "subtypes011", oldID = "DS_Epithelial_Cells_1", newID = "Cholangiocytes",
               newOrder = gsub("DS_Epithelial_Cells_1","Cholangiocytes",levels(Mobj$subtypes011)))
Mobj=modifyAno(Mobj, anoCol = "subtypes011", oldID = "RPS_Epithelial_Cells_2", newID = "Granulosa",
               newOrder = gsub("RPS_Epithelial_Cells_2","Granulosa",levels(Mobj$subtypes011)))
Mobj=modifyAno(Mobj, anoCol = "subtypes011", oldID = "RPS_CS_Epithelial_Cells_1", newID = "Epi_undef1",
               newOrder = gsub("RPS_CS_Epithelial_Cells_1","Epi_undef1",levels(Mobj$subtypes011)))
Mobj=modifyAno(Mobj, anoCol = "subtypes011", oldID = "MIX_Epithelial_Cells_1", newID = "Epi_undef2",
               newOrder = gsub("MIX_Epithelial_Cells_1","Epi_undef2",levels(Mobj$subtypes011)))
Mobj=modifyAno(Mobj, anoCol = "subtypes011", oldID = "Xendo_unk.", newID = "Mast",
               newOrder = gsub("Xendo_unk.","Mast",levels(Mobj$subtypes011)))
Mobj$subtypes011=factor(Mobj$subtypes011, levels = c(
  levels(Mobj$subtypes011)[1:30],"Mast",levels(Mobj$subtypes011)[31:71],levels(Mobj$subtypes011)[73:107]
))

Mobj_tmp=subset(Mobj, subtypes011 %!in% c("Mast","Prolif."))
Mobj_tmp$subtypes011=droplevels(Mobj_tmp$subtypes011)
xxx=table(Mobj_tmp$System, Mobj_tmp$subtypes011)
yyy=colPercents(xxx)[1:dim(xxx)[1],]
yyy=melt(yyy)

DEG_clusterUniq=function(tmpmarkers, minGnum=10){

final.clusters <- tmpmarkers %>% dplyr::count(cluster) %>% filter(n>minGnum) %>% pull(cluster)
markers.final <- tmpmarkers %>% filter(cluster %in% final.clusters)
# Remove duplicated genes
markers.final$rows <- 1:nrow(markers.final)
dup <- markers.final[duplicated(markers.final$gene), ]$gene
remove <- purrr::map(.x=dup, .f=function(i){
  keep <- markers.final %>%
    dplyr::filter(gene== i) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    head(1) %>%
    pull(rows)

  remove <- markers.final %>%
    filter(gene==i) %>%
    filter(rows!=keep) %>%
    pull(rows)

  return(remove)
}) %>% unlist()
if(!is.null(remove)){
markers.final <- markers.final[-c(remove), ]
}else{
  markers.final <- markers.final
}
final.clusters <- markers.final %>% dplyr::count(cluster) %>% filter(n>minGnum) %>% pull(cluster)
markers.final <- markers.final %>% filter(cluster %in% final.clusters)
rownames(markers.final) <- markers.final$gene

return(markers.final)
}

topDEGrun=function(DEGs, topn=10, clus=NULL){
  DEGs$pct.diff=DEGs$pct.1-DEGs$pct.2
  DEGs=as.data.frame(DEGs %>% group_by(cluster) %>% top_n(n=topn, wt=(avg_log2FC*pct.diff)))
  if(is.null(clus)){
    clus=unique(DEGs$cluster)
  }
  return(DEGs[DEGs$cluster %in% clus,])
}

AllDEG2CorPlot=function(obj,AllDEG,identX="seurat_clusters",topn=50,filter=FALSE,rmv=NULL,log_weight=1){

  DevNeuDEGold_top50=as.data.frame(AllDEG %>% group_by(cluster) %>% top_n(n=topn, wt=avg_log2FC*(log_weight+(1-log_weight)*(pct.1-pct.2))))
  xxx=DevNeuDEGold_top50
  if(filter){
    #xxx=xxx[(!grepl("^IG[LHK]", xxx$gene)),]
    #xxx=xxx[(!grepl("^HIST", xxx$gene)),]
    xxx=xxx[(!grepl("^RP[LS]", xxx$gene)),]
    xxx=xxx[(!grepl("^MT-", xxx$gene)),]
  }
  DevNeuDEGold_top50=xxx
  NeuCorData=obj
  Idents(NeuCorData)=identX
  tdata1 <- AverageExpression(subset(NeuCorData, seurat_clusters %!in% rmv),
                              features = unique(DevNeuDEGold_top50$gene),
                              return.seurat = TRUE)
  Mat_Tmp <- as.matrix(tdata1@assays$RNA$data)
  ComplexHeatmap::pheatmap(cor(Mat_Tmp),  name ="correlation", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2", clustering_distance_rows = "correlation",
                           clustering_distance_cols = "correlation") #clustering_method = "complete", clustering_distance_rows = "euclidean", #, treeheight_col = 0, treeheight_row = 0)  #scale = "row",

}

DEGs_Mobj_107subCT=data.frame()
Idents(Mobj_keep)="subtypes015"
for(loc in levels(Mobj_keep$Location_short)){
  Mobj_keep_tmp=subset(Mobj_keep, Location_short %in% c(loc))
  Mobj_keep_tmp$Location_short=droplevels(Mobj_keep_tmp$Location_short)
  Mobj_keep_tmp$subtypes015=droplevels(Mobj_keep_tmp$subtypes015)

  myCT=names(table(Mobj_keep_tmp$subtypes015)[table(Mobj_keep_tmp$subtypes015)>=5])
  if(length(myCT)>=1){
  for(x in myCT){
    tmp_107subCT=FindMarkers(Mobj_keep_tmp, ident.1 = x, group.by = "subtypes015", logfc.threshold = 0.25, min.pct = 0.2, only.pos = T, max.cells.per.ident = 2000, test.use = "t")
    if(dim(tmp_107subCT)[1]>=3){
    tmp_107subCT$gene=rownames(tmp_107subCT)
    tmp_107subCT$subtypes=x
    tmp_107subCT$location=loc

    DEGs_Mobj_107subCT=rbind(DEGs_Mobj_107subCT, tmp_107subCT)
    }
  }
  }
}
DEGs_Mobj_107subCT=DEGs_Mobj_107subCT[DEGs_Mobj_107subCT$p_val_adj<=0.05,]
