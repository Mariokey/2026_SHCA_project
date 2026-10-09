# ============================================================================
# Gene co-expression modules with hdWGCNA
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Identification of cellular modules"
# Source     : extracted from Analysis240525.r (L5846-6222, L6892-7214)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
subdata=FindVariableFeatures(subdata, method = "vst", nfeatures=2000)
subdata=ScaleData(subdata, assay = "RNA",
                  vars.to.regress = "pct_counts_mt",
                  verbose = FALSE)
cmplist=NULL #list(c("Artery","Spleen"))

wgcna_MP=wgcnaRun(subdata, baseDir="./hdWGCNA/", harmony="pca_harmony",
                  prefix = "MP", fraction=0.1,
                  submode=TRUE, cmplist=cmplist,
                  ctname="MYL_SubClass02",
                  ct_vec = names(table(subdata$MYL_SubClass02)[table(subdata$MYL_SubClass02)>120]),
                  gpname="Organ")

sc=readRDS("./hdWGCNA/MP.hdWGCNA_object.rds")
wgcna_MP=wgcnaRun(sc, baseDir="./hdWGCNA/", harmony="pca_harmony",
                  prefix = "MP", fraction=0.1,
                  submode=TRUE, cmplist=cmplist,
                  ctname="MYL_SubClass02",
                  ct_vec = names(table(subdata$MYL_SubClass02)[table(subdata$MYL_SubClass02)>120]),
                  gpname="Organ", gp_vec = names(table(sc$Organ)[table(sc$Organ)>0]),
                  plot_mode = TRUE)

png(paste0("./hdWGCNA/","MP","_hdWGCNA_Lollipop.png"), res = 300, width = 21*300, height = 21*300)
plot_grid(plotlist = wgcna_MP$ModuleRadarPlot, ncol = 3, align = "hv")
dev.off()

png(paste0("./hdWGCNA/","MP","_ModuleDotPlot_split.png"), res = 300, width = 24*300, height = 6*300)
wgcna_MP$ModuleDotPlot_split
dev.off()

write.xlsx(wgcna_MP$MEfun_GO, file = paste0("./hdWGCNA/","MP","_ModuleGOfunPlot.xlsx"))

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

plot2dfun=function(obj,features=c("nFeature_RNA","percent.mt","MKPP_score1")){

  plotlist=list()
  d2=features[1:2]
  dothers=features[3:length(features)]
  for(x in dothers){
    metaplotDF=obj@meta.data
    if(x %!in% colnames(metaplotDF)){
      tmpvplot=VlnPlot(obj, features = x)
      metaplotDF=cbind(metaplotDF,tmpvplot[[1]]$data[rownames(metaplotDF),])
    }
    xxx=ggplot(metaplotDF[order(metaplotDF[, x], decreasing = FALSE),])+
      geom_point(aes_string(x= d2[1], y= d2[2], fill=x), shape=21, size=1)+
      scale_fill_continuous(type = "viridis")+scale_color_continuous(type = "viridis")+
      theme_cowplot()

    plotlist[[x]]=xxx
  }

  return(plotlist)
}

runHeat<-function(t_data, qt=c(0.9,0.95), targetCT="HSC", type=FALSE, ...){

tmpdf=t_data@meta.data[,c("Location_short","Mono_score.01","MoDM_score.01","FeDM_score.01")]
tmpdf[,2]=rescale(tmpdf[,2], to = c(0,1))
tmpdf[,3]=rescale(tmpdf[,3], to = c(0,1))
tmpdf[,4]=rescale(tmpdf[,4], to = c(0,1))
AA=tmpdf[,c(1,2)] %>% group_by(Location_short) %>% summarise_at("Mono_score.01",mean) %>% as.data.frame()
BB=tmpdf[,c(1,3)] %>% group_by(Location_short) %>% summarise_at("MoDM_score.01",mean) %>% as.data.frame()
CC=tmpdf[,c(1,4)] %>% group_by(Location_short) %>% summarise_at("FeDM_score.01",mean) %>% as.data.frame()
AABBCC=left_join(AA, left_join(BB,CC, by="Location_short"), by="Location_short")
rownames(AABBCC)=AABBCC$Location_short
ggplot(AABBCC)+geom_point(aes(Mono_score.01, FeDM_score.01, color=Location_short))+
  geom_text(data = AABBCC, aes(Mono_score.01,FeDM_score.01,label=Location_short, color=Location_short), nudge_x = 0.01, nudge_y=0.01, check_overlap = T, size=4, family="Helvetica")+
  theme_cowplot()+NoLegend()
p=ggplot(AABBCC)+geom_point(aes(MoDM_score.01, FeDM_score.01, color=Location_short))+
  geom_text(data = AABBCC, aes(MoDM_score.01,FeDM_score.01,label=Location_short, color=Location_short), nudge_x = 0.01, nudge_y=0.01, check_overlap = T, size=4, family="Helvetica")+
  theme_cowplot()+NoLegend()
ggsave(paste0("./",'SX_MP_f.',"cor_1",'.pdf'),p,width = 5.5, height = 4)
p=ggplot(AABBCC)+geom_point(aes(MoDM_score.01, Mono_score.01, color=Location_short))+
  geom_text(data = AABBCC, aes(MoDM_score.01,Mono_score.01,label=Location_short, color=Location_short), nudge_x = 0.01, nudge_y=0.01, check_overlap = T, size=4, family="Helvetica")+
  theme_cowplot()+NoLegend()
ggsave(paste0("./",'SX_MP_f.',"cor_2",'.pdf'),p,width = 5.5, height = 4)

Mono_score=c("FCN1","VCAN","S100A12","FYN","CD300E","CCR2","CLU","CRISPLD2")  #"LILRB1",CD14","FCGR3A","MARC1","CRISPLD2"
MoDM_score=c("CCL20","C1QA","APOE","PHLDA3","FABP4","CCL18","CD52","BCL22A1A","S100A4","S100A6","S100A11","S100A13","CORO1A","FXDY5","HLA-DRA","HLA-DRB1","S100A8","S100A9","TIMP1","MARCO","TMSB10","HLA-DPB1","HLA-DPA1","HLA-DQA2")
FeDM_score=c("FOLR2","TIMD4","RGS1","STARD13","CD163","LYVE1","NINJ1","GAS6","MRC1","CD163","MERTK","ABCA1","SIGLEC1","CSF1R","SELENOP","PTLP")

t_data=AddModuleScore(t_data, features = list(Mono_score[1:5]),
                     ctrl = 500, name = "Mono_score.0")
t_data=AddModuleScore(t_data, features = list(MoDM_score),
                     ctrl = 500, name = "MoDM_score.0")
t_data=AddModuleScore(t_data, features = list(FeDM_score),
                     ctrl = 500, name = "FeDM_score.0")

  qt=c(0.35,1)
  type=TRUE

  targetCT=c("Mono","MoDM","FeDM")

Mono_score=c("VCAN","FCN1","S100A12","S100A8")
MoDM_score=c("TIMP1","S100A4","S100A6","S100A11","MARCO") ###"HLA-DPB1","HLA-DPA1"#c("CD52","BCL22A1A","S100A4","S100A6","S100A11","FXDY5","HLA-DRA","HLA-DRB1","S100A8","S100A9","TIMP1","MARCO","TMSB10","HLA-DPB1","HLA-DPA1","HLA-DQA2")
FeDM_score=c("PLTP","FOLR2","LYVE1","STARD13","CD163","SELENOP")  ###,"CD163" "TIMD4","RGS1",
STEMsigG=Mono_score
CMPsigG=MoDM_score
GMPsigGenes=FeDM_score
LYMsig1=c()

  heatmapGenes=unique(c(STEMsigG,
                        CMPsigG,
                        GMPsigGenes,
                        LYMsig1))

gA=c()
heatmapGenesx=c()
for(i in 1:length(heatmapGenes)){
  g=heatmapGenes[i]
  if(g %in% rownames(t_data)){
    gA=c(gA, ifelse(g %in% STEMsigG, "Mono",
                    ifelse(g %in% CMPsigG, "MoDM",
                          ifelse(g %in% GMPsigGenes, "FeDM", "Udf"))))
    heatmapGenesx=c(heatmapGenesx,g)
   }
}
heatmapGenes=heatmapGenesx

myanno=data.frame(gene=heatmapGenes, CellSign=gA)
myanno1=data.frame(CellSign=myanno$CellSign)
myanno$CellSign=factor(myanno$CellSign, levels = c("Mono","MoDM","FeDM","Udf"))
rownames(myanno1)=myanno$gene

myanno1$ccc=rownames(myanno1)
myanno1=myanno1[myanno1$CellSign %in% targetCT,]
myanno1$ccc=NULL

ann_colors = list(
  CellSign=c(Mono = "blue", MoDM = "#E7298A", FeDM = "#66A61E")## CLP = "red")
)

  heatmapGenes=rownames(myanno1)

t_data_P=t_data
t_data_H=t_data

qL=quantile(t_data_P$Mono_score.01, probs = qt[1])
qH=quantile(t_data_P$Mono_score.01, probs = qt[2])
t_data_P_qt=subset(t_data_P, (Mono_score.01>=qL[[1]] & Mono_score.01<=qH[[1]]))
qL=quantile(t_data_H$FeDM_score.01, probs = qt[1])
qH=quantile(t_data_H$FeDM_score.01, probs = qt[2])
t_data_H_qt=subset(t_data_H, (FeDM_score.01>=qL[[1]] & FeDM_score.01<=qH[[1]]))

  myColor <- colorRampPalette(c(rep("darkblue",1), rep("white", 1), rep("#A50000", 1)))(100)

heatmapGenes00=heatmapGenes

hmpH <- subset(t_data_H_qt,features = heatmapGenes00)
Mat_Tmp <- as.matrix(hmpH@assays$RNA$data)
colnames(Mat_Tmp) <- paste0(colnames(Mat_Tmp))
Xmtx1=cor(t(Mat_Tmp), method = "spearman")   ###"pearson", "kendall", "spearman"
Xmtx1=Xmtx1[heatmapGenes00, heatmapGenes00]
range <- max(abs(Xmtx1), na.rm = TRUE)
tbreak = seq(-range/1.5, range/1.5, length.out = 150)

HHH=ComplexHeatmap::pheatmap(Xmtx1, legend = F, na_col = "grey90", color = myColor, breaks = tbreak, fontsize = 18, #4.5
                             cluster_rows = F, cluster_cols = F, #treeheight_col = 0, treeheight_row = 0,
                             border_color = "white", name = "Cor.", show_rownames = T, show_colnames = F,
                             right_annotation  = rowAnnotation(df = data.frame(CellSign=myanno1[rownames(Xmtx1),], row.names = rownames(Xmtx1)), col=ann_colors, show_legend=F, show_annotation_name = F),
                             row_names_side = "left") #, annotation_colors = ann_colors)
Horder=HHH@row_names_param$labels[row_order(draw(HHH))]

hmpP <- subset(t_data_P_qt,features = heatmapGenes00)
Mat_Tmp <- as.matrix(hmpP@assays$RNA$data)
colnames(Mat_Tmp) <- paste0(colnames(Mat_Tmp))
Xmtx=cor(t(Mat_Tmp), method = "spearman")
Xmtx=Xmtx[heatmapGenes00, heatmapGenes00]
range <- max(abs(Xmtx), na.rm = TRUE)
tbreak = seq(-range/1.5, range/1.5, length.out = 150)

Porder=Horder[Horder %in% rownames(Xmtx)]
addOrder=rownames(Xmtx)[rownames(Xmtx) %!in% Horder]
Porder=c(addOrder, Porder)
PPP=ComplexHeatmap::pheatmap(Xmtx[Porder,Porder], legend = T, na_col = "grey90", color = myColor, breaks = tbreak, fontsize = 30, #4.5
                             cluster_rows = F, cluster_cols = F, #treeheight_col = 0, treeheight_row = 0,
                             border_color = "white", name = "Cor.", show_rownames = F, show_colnames = F,
                             right_annotation  = rowAnnotation(df = data.frame(CellSign=myanno1[Porder,], row.names = Porder), col=ann_colors, show_annotation_name = F),
                             row_names_side = "left") #, annotation_colors = ann_colors)

  HP=plot_grid(HHH, PPP)
  t_file=paste0(out0203,"Fig.HeatmapGeneModule.", paste(qt[1], qt[2], sep = "_"),".HSC_HC_MS00")
  runFig(t_file, HP, 32,18)
  addPPT(my_ppt, paste0("Gene Modules enriched from HSC populations (left: Healthy; right: Patients) with specified stemness score range of ", paste(qt[1], qt[2], sep = " to ")), t_file, type=TRUE)

  #  }
}

wgcnaRun=function(sc, prefix="tmp", fraction=0.1, submode=TRUE, cmplist=NULL, harmony=NULL, min_cells = 50,
                  ctname="seurat_clusters", ct_vec=NULL, gpname="Sample", gp_vec=NULL, baseDir="./", plot_mode=FALSE){

  library(devtools)
  library(GeneOverlap)
  library(UCell)
  library(ggrepel)
  library(igraph)
  library(ggforestplot)
  library(WGCNA)
  library(hdWGCNA)
  library(Seurat)
  library(cowplot)
  library(patchwork)

  dir.create(baseDir, recursive = T)

  set.seed(12345)
  theme_set(theme_cowplot())

  ### subtypes01 majortypes Sample

  GmodRes=list()

  harmonized=FALSE
  if(!is.null(harmony)){
    harmonized=TRUE
  }

  if(!plot_mode){

  if(!is.null(ct_vec)){
    sc=sc[,sc@meta.data[,ctname] %in% ct_vec]
    sc@meta.data[,ctname]=droplevels(sc@meta.data[,ctname])
  }
  if(!is.null(gp_vec)){
    sc=sc[,sc@meta.data[,gpname] %in% gp_vec]
    sc@meta.data[,gpname]=droplevels(sc@meta.data[,gpname])
  }

  print(table(sc@meta.data[,ctname],sc@meta.data[,gpname]))

  if(fraction>0){
    sc <- SetupForWGCNA(
      sc,
      gene_select = "fraction",
      fraction = fraction, # fraction of cells that a gene needs to be expressed in order to be included
      group.by=ctname,
      wgcna_name = paste0(prefix,"_hdWGCNA") # the name of the hdWGCNA experiment
    )
  }else{
    if(submode){
      sc=FindVariableFeatures(sc, method = "vst", nfeatures=2000)
    }
    sc <- SetupForWGCNA(
      sc,
      gene_select = "variable", # the gene selection approach
      wgcna_name = paste0(prefix,"_hdWGCNA") # the name of the hdWGCNA experiment
    )
  }
  reduc_name="pca"
  if(!is.null(harmony)){
    reduc_name=harmony
  }
  sc <- MetacellsByGroups(
    seurat_obj = sc,
    group.by = c(ctname, gpname), # specify the columns in seurat_obj@meta.data to group by
    reduction = reduc_name, # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    min_cells = min_cells,
    ident.group = ctname # set the Idents of the metacell seurat object
  )
  sc <- NormalizeMetacells(sc)
  #####################################################

  ct_vec=table(sc@misc[[paste0(prefix,"_hdWGCNA")]]$wgcna_metacell_obj[[ctname]][[1]]) %>% names()
  gp_vec=table(sc@misc[[paste0(prefix,"_hdWGCNA")]]$wgcna_metacell_obj[[gpname]][[1]]) %>% names()
  sc <- SetDatExpr(
    sc,
    group_name = ct_vec, # the name of the group of interest in the group.by column
    group.by=ctname, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data', # using normalized data
    multi.group.by = gpname,
    multi_group_name = gp_vec
  )

  sc <- TestSoftPowers(
    sc,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )

  sc <- ConstructNetwork(
    sc,
    overwrite_tom = TRUE,
    tom_outdir = baseDir,
    tom_name = prefix # name of the topoligical overlap matrix written to disk
  )

  # need to run ScaleData first or else harmony throws an error:

  # compute all MEs in the full single-cell dataset
  harmony_var=NULL
  if(!is.null(harmony)){
    harmony_var="Sample"
  }
  sc <- ModuleEigengenes(
    sc,
    group.by.vars=harmony_var # No harmony or #"Sample"
  )

  sc <- ModuleConnectivity(
    sc,
    group.by = ctname, harmonized = harmonized, group_name = ct_vec,
  )
  sc <- ResetModuleNames(
    sc,
    new_name = paste0(prefix,"-M")
  )

  #compute gene scoring for the top 25 hub genes by kME for each module
  sc <- ModuleExprScore(
    sc,
    n_genes = 25,
    method='UCell'
  )
  # get hMEs from seurat object
  t_MEs <- GetMEs(sc, harmonized=harmonized)
  # add hMEs to Seurat meta-data:
  sc@meta.data <- cbind(sc@meta.data, t_MEs)

  saveRDS(sc, file=paste0(baseDir,prefix,'.hdWGCNA_object.rds'))

  }
  ################################################################################

tmpTOM <- GetTOM(sc)
GmodRes[[paste0("TOM_matrix")]]=tmpTOM

MEs <- GetMEs(sc, harmonized=harmonized) #hMEs <- GetMEs(seurat_obj)
MEs_names=sort(colnames(MEs)); print(MEs_names)
GmodRes[[paste0("MEs_df")]]=MEs

modules <- GetModules(sc) %>% subset(module != 'grey') ###greys SHOULD be removed in ALL downstream analysis
# show the first 6 columns:
head(modules)
GmodRes[[paste0("modules_df")]]=modules

# get top25 hub genes
hub_df <- GetHubGenes(sc, n_hubs = 100)
head(hub_df)
GmodRes[[paste0("hub100genes_df")]]=hub_df

hub_df$gene=hub_df$gene_name
hub_df$avg_log2fc=100
hub_df$cluster=hub_df$module
ResFun=enrichFUNgo(ud_table = hub_df[hub_df$cluster %!in% c("grey"),],sp="hsa",intermode = FALSE, enriches = c("GO"))
GmodRes[[paste0("MEfun_GO")]]=ResFun$dflist$GO
GOdf=ResFun$dflist$GO
ddd=GOdf%>% group_by(celltype) %>% top_n(n=5, wt=-pvalue)
ddd$celltype=factor(ddd$celltype, levels = MEs_names[MEs_names %!in% c("grey")])
ddd$Description=factor(ddd$Description, levels = rev(unique(ddd$Description)))

GmodRes[[paste0("MEfun_GOplot")]]=ggplot(ddd)+
  geom_point(aes(celltype, Description, color=-log10(pvalue), size=Factor))+
  xlab("")+ylab("")+theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(paste0(baseDir,prefix,"_hdWGCNA_Dendrogram.pdf"), width = 8, height = 4)
PlotDendrogram(sc, main=paste0(prefix,' hdWGCNA Dendrogram'))
dev.off()

GmodRes[["KMEsPlot"]]=PlotKMEs(sc, ncol=4)

GmodRes[["ModuleRadarPlot"]]=ModuleRadarPlot(
  sc[,sc@meta.data[, ctname] %in% ct_vec],
  group.by = ctname,
  barcodes = sc[,sc@meta.data[, ctname] %in% ct_vec] %>% colnames(),
  axis.label.size=4,
  grid.label.size=4,
  combine=FALSE
)

## make a featureplot of hMEs for each module
#)
## stitch together with patchwork

# plot with Seurat's DotPlot function
t_modules <- GetModules(sc)
t_mods <- levels(t_modules$module); t_mods <- t_mods[t_mods != 'grey']
#xxx=DotPlot(sc, features=t_mods, group.by = ctname)+
GmodRes[["ModuleDotPlot_split"]]=splitDotPlot(sc[,sc@meta.data[, ctname] %in% ct_vec],
                                              group.by=ctname, groupV=NULL, split.by=gpname, splitV=gp_vec, featureV=t_mods)

  if(!is.null(cmplist)){

    # list of clusters to loop through
    clusters <- ct_vec
    cmplist=cmplist #list(c("B-A-3D","B"))
    # set up an empty dataframe for the DMEs
    DMEs <- data.frame()
    # all clusters as one
    for(i in 1:length(cmplist)){
      cmpitem=paste(cmplist[[i]], collapse = " vs. ", sep="")

      # identify barcodes for group1 and group2 in eadh cluster
      group1 <- sc[,sc@meta.data[, gpname] %in% cmplist[[i]][1] ] %>% colnames()
      group2 <- sc[,sc@meta.data[, gpname] %in% cmplist[[i]][2] ] %>% colnames()

      # run the DME test
      cur_DMEs <- FindDMEs(
        sc,
        barcodes1 = group1,
        barcodes2 = group2,
        test.use='wilcox',
        pseudocount.use=0.01, # we can also change the pseudocount with this param
        wgcna_name = sc@misc$active_wgcna #'INH'
      )

      # add the cluster info to the table
      cur_DMEs$cluster <- "All"
      cur_DMEs$cmp <- cmpitem

      # append the table
      DMEs <- rbind(DMEs, cur_DMEs)

    }
    # loop through the clusters
    for(cur_cluster in clusters){

      for(i in 1:length(cmplist)){
        cmpitem=paste(cmplist[[i]], collapse = " vs. ", sep="")

        # identify barcodes for group1 and group2 in eadh cluster
        group1 <- sc[,sc@meta.data[, ctname] %in% cur_cluster &
                       sc@meta.data[, gpname] %in% cmplist[[i]][1] ] %>% colnames()
        group2 <- sc[,sc@meta.data[, ctname] %in% cur_cluster &
                       sc@meta.data[, gpname] %in% cmplist[[i]][2] ] %>% colnames()

        if(length(group1)<5 | length(group2)<5) next

        # run the DME test
        cur_DMEs <- FindDMEs(
          sc,
          barcodes1 = group1,
          barcodes2 = group2,
          test.use='wilcox',
          pseudocount.use=0.01, # we can also change the pseudocount with this param
          wgcna_name = sc@misc$active_wgcna #'INH'
        )

        # add the cluster info to the table
        cur_DMEs$cluster <- cur_cluster
        cur_DMEs$cmp <- cmpitem

        # append the table
        DMEs <- rbind(DMEs, cur_DMEs)
      }
    }

    GmodRes[["InterPlot_split"]]=PlotDMEsLollipop(
      sc, #[,sc@meta.data[,ctname] %in% unique(DMEs$cluster)],
      DMEs,
      wgcna_name=sc@misc$active_wgcna,
      group.by = "cluster",
      comparison = c("All",ct_vec[ct_vec %in% unique(DMEs$cluster)]),
      pvalue = "p_val_adj"
    )

  }

  return(GmodRes)

}
