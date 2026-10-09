# ============================================================================
# Public dataset cell-type projection with scibet and CNTN4+ epithelial cell recovery
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Public scRNA-seq data processing"
# Source     : extracted from tmp20250610.r (L1-815)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
source("DATA_DIR")
setwd("./data/")

pub_data_list=readRDS("./data/PubData/Epithelial_list.rds")
print(length(pub_data_list))
pub_data=merge(pub_data_list[[1]],
               pub_data_list[c(2,4:7,9,
                               11:35,
                               37:38,
                               41:length(pub_data_list))])
saveRDS(pub_data, paste0("./PubData/","Epithelial",".rds"))

data_genes=read.csv("./PubData/Epithelial_feature.csv")
pcgenes=read.csv("./ref/hum_pcgene.list", header = F)
mygenes=intersect(pcgenes$V1, data_genes$features)

pub_data=subset(pub_data, features = mygenes)
pub_data=UpdateSeuratObject(pub_data)
pub_data=subset(pub_data, nFeature_RNA>=200)
saveRDS(pub_data, "./data/PubData/Epithelial_final.rds")

###
xxx001=readRDS("DATA_DIR")
xxx002=readRDS("DATA_DIR") ### #1
#meta_CellOClassification meta_phase    meta_scibetHCL           meta_scibetHPA meta_scibetMajorType meta_scibetSubType
subpub_data=subset(pub_data, cellID %in% colnames(xxx002))
pheatmap::pheatmap(log2(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)+1))
tmpmat=log2(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)+1)

cn02=c("Enterocytes","BEST4_Epithelial_Cells","Pit_Cells","Chief_Cells","Papillary_Epithelial_Cells","Mid_Suprabasal_Cells","Late_Suprabasal_Cells","Multiciliated_Cells","Hepatocytes","Goblet_Cells","Endocrine_Cells")
rn02=c("Intestinal epithelium",
"Gastrointestinal epithelium",
"Distal lung epithelium",
"Basal like",
"Ciliated",
"Hepatocyte",
"MUC2+ goblet",
"Enteroendocrine")
mat02=colPercents(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)[rn02,cn02])[1:length(rn02),]
ppp=pheatmap::pheatmap(mat02[rn02,cn02], cluster_rows = F, cluster_cols = F, color = viridis::viridis_pal()(100))
ggsave(filename = "./data/NO2_pub_match.pdf", plot = ppp,  width = 5.7, height = 4.1)

cn01=c("Spinous_Cells","Early_Suprabasal_Cells","Basal_Cells","Vellus_Hair_Follicle","Cholangiocytes","Enterocytes","Club_Cells","Pit_Cells","Tuft_Cells","Chief_Cells","Colonocytes","BEST4_Epithelial_Cells","Mid_Suprabasal_Cells","Airway_Suprabasal_Cells")
rn01=c("Spinous Epithelial Cell KRT1_high","Keratinocyte KRT1_high","Basal Epithelial Cell POSTN_high","Follicular Epithelial Cell","Cholangiocyte HIST1H2AM_high","Enterocyte APOA1_high","Secretory Cell","Pit Mucosal Epithelial Cell","Tuft Cell","Epithelial Cell TFF3_high","Simple Epithelial Cell","Absorptive Cell","Squamous Epithelial Cell KRT13_high","Basal Cell KRT17_high")
setdiff(rn01,rownames(tmpmat))
setdiff(cn01,colnames(tmpmat))
mat01=colPercents(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)[rn01,cn01])[1:length(rn01),]
ppp=pheatmap::pheatmap(mat01[rn01,cn01], cluster_rows = F, cluster_cols = F, color = viridis::viridis_pal()(100))
ggsave(filename = "./data/NO1_pub_match.pdf", plot = ppp,  width = 8, height = 6)

##################################################################################################################################################

pub_data=readRDS(paste0("./PubData/","Immune",".rds"))
data_genes=read.csv("./PubData/Immune_feature.csv")
mygenes=intersect(pcgenes$V1, data_genes$features)

pub_data=subset(pub_data, features = mygenes)
pub_data=UpdateSeuratObject(pub_data)
pub_data=subset(pub_data, nFeature_RNA>=200)
saveRDS(pub_data, "./data/PubData/Immune_final.rds")
dim(pub_data)

pub_data_imm=readRDS("./workspace/PRJ21_mOrganYang/R_data/PubData/Immune_final_pred.rds")
xxx001=readRDS("DATA_DIR") ### #1
subpub_data=subset(pub_data_imm, cellID %in% colnames(xxx001))
pheatmap::pheatmap(log2(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)+1))
tmpmat=log2(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)+1)

cn02=colnames(tmpmat)[c(1:7,20:22,13:19)]
rn02=c("B cell",
       "T cell/NK cell 2",
       "Macrophage/monocyte 2")
mat02=colPercents(table(subpub_data@meta.data$clusterName,subpub_data@meta.data$label)[rn02,cn02])[1:length(rn02),]
ppp=pheatmap::pheatmap(mat02[rn02,cn02], cluster_rows = F, cluster_cols = F, color = viridis::viridis_pal()(100))
ggsave(filename = "./data/NO2_pub_match_IMM.pdf", plot = ppp,  width = 8.8, height = 3.5)

#

#########################
pub_data=readRDS("./PubData/Epithelial_final_pred.rds")
CNTN4_pubdata=subset(pub_data,label %in% c("CNTN4_NEL_Cells"))
sample_gt30=names(table(CNTN4_pubdata$sampleID)[table(CNTN4_pubdata$sampleID)>=10])   ### 117 samples and 11709 cells
CNTN4_pubdata=subset(CNTN4_pubdata, sampleID %in% c(sample_gt30))
CNTN4_pubdata=reDRClusterTest(CNTN4_pubdata, sp="hum", ndim=20, nClu=100, mdist=0.25,
                formalFeature=NULL, mypca="harmony", harmony=10, prepro=TRUE,
                mybatch="sampleID", umap=TRUE, vstN=2000, regVar=NULL, rmvHVGribo=FALSE)  ###c("stats_mitoPercent")

DimPlot(CNTN4_pubdata, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(CNTN4_pubdata, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
scPioneer::DimPlot_idx(CNTN4_pubdata,  group.by = "seurat_clusters", reduction = "umap")

CNTN4_data$System=factor(CNTN4_data$System, levels = levels(Mobj_keep$System))
CNTN4_data$Location_short=factor(CNTN4_data$Location_short, levels = levels(Mobj_keep$Location_short))
inobj=CNTN4_data
inobj$UMAP_1=inobj@reductions$umap@cell.embeddings[,1]
inobj$UMAP_2=inobj@reductions$umap@cell.embeddings[,2]
gpColor=hue_pal()(length(levels(inobj$System)))
names(gpColor)=levels(inobj$System)
plotlist=list()
for(x in names(gpColor)){
  ppp=ggplot()+stat_density_2d(data=inobj@meta.data, #%>% sample_n(size=10000),
                               aes(UMAP_1, UMAP_2, color=as.factor("..level..")), breaks=c(0.0005,100), h = c(0.8, 0.8), n=100)+
    scale_color_manual(values = c("grey",rep("NA",100),"NA"))+
    stat_density_2d(data=inobj@meta.data[inobj@meta.data$System %in% x, ],
                    aes_string(x = "UMAP_1", y = "UMAP_2",  alpha = "..level..", fill = "System"),
                    linewidth = 0, geom = "polygon",  n = 200, h = c(1, 1)) + theme_cowplot() +
    scale_fill_manual(values = gpColor)+
    scale_alpha_continuous(limits=c(0,0.3)) + NoLegend()
  plotlist[[x]]=ppp
}
plot_grid(plotlist = plotlist, ncol = 3, labels = names(plotlist))

############## CNTN4>0
CNTN4_pubdata2=subset(pub_data,CNTN4>0)
sample_gt30_2=names(table(CNTN4_pubdata2$sampleID)[table(CNTN4_pubdata2$sampleID)>=10])   ### 117 samples and 11709 cells
CNTN4_pubdata2=subset(CNTN4_pubdata2, sampleID %in% c(sample_gt30_2))
CNTN4_pubdata2=reDRClusterTest(CNTN4_pubdata2, sp="hum", ndim=20, nClu=100, mdist=0.25,
                              formalFeature=NULL, mypca="harmony", harmony=10, prepro=TRUE,
                              mybatch="sampleID", umap=TRUE, vstN=2000, regVar=NULL, rmvHVGribo=FALSE)  ###c("stats_mitoPercent")

DimPlot(CNTN4_pubdata2, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(CNTN4_pubdata2, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
scPioneer::DotPlot2(CNTN4_pubdata2, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"))

hhh=intersect(colnames(CNTN4_pubdata2), colnames(subset(CNTN4_pubdata, seurat_clusters %in% c(2,5,7,8,9,10))))
DimPlot(CNTN4_pubdata2, cells.highlight = hhh)

Marker_Lev2=c(
  "MS4A1","MZB1",
  "CD3D","KLRD1",
  "CD1C","FCN1","C1QC","FCGR3B","CPA3",

  "PGC","PHGR1","PIGR", #Gastric/Intestinal
  "NEUROD1","RFX6",  #Endocrine
  "PRSS1","PRSS2",  #Pancreas
  "ALB","HP","KRT19","CD24", #Hep/BEC  #"CFTR","SOX9"
  "CHAT","POU2F3", #Tuft
  "KRT5","KRT15","KRT1","KRT10", #Basal/Super
  "CALML3", #Skin
  "WFDC2", #Tracheal
  "KRT7","CLDN4", #Urinary_tract
  "MAL","CDH16", #Collecting-ductal(CCDPC/Kidney)
  "SERPINE2","PRKAR2B", #Granulosa #GSTA1
  "PRG4","WT1", #Mesothelial
  "CNTN4", "PDZRN3", #NE_like "YAP1","GPC6", "FBXL7",

  "EGFL7","COL1A2","ACTA2",
  "TCF3","SYCP3","TNP1",  #TEX15/MAGEB2
  "MKI67"
)
scPioneer::DotPlot2(CNTN4_pubdata2,features = Marker_Lev2, group.by = "seurat_clusters")

hhh2=c(colnames(CNTN4_pubdata2), colnames(CNTN4_pubdata)) %>% unique()
CNTN4_pubdata3=subset(pub_data, cb %in% hhh2)

CNTN4_pubdata3=reDRClusterTest(CNTN4_pubdata3, sp="hum", ndim=20, nClu=100, mdist=0.25,
                               formalFeature=NULL, mypca="harmony", harmony=10, prepro=TRUE,
                               mybatch="sampleID", umap=TRUE, vstN=2000, regVar=NULL, rmvHVGribo=FALSE)  ###c("stats_mitoPercent")

DimPlot(CNTN4_pubdata3, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(CNTN4_pubdata3, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
scPioneer::DotPlot2(CNTN4_pubdata3, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"))
scPioneer::DotPlot2(CNTN4_pubdata3,features = Marker_Lev2, group.by = "seurat_clusters")
sample_gt30_3=names(table(CNTN4_pubdata3$sampleID)[table(CNTN4_pubdata3$sampleID)>=10])   ### 117 samples and 11709 cells
CNTN4_pubdata3=subset(CNTN4_pubdata3, sampleID %in% c(sample_gt30_3))
CNTN4_pubdata3_sub=subset(CNTN4_pubdata3, meta_tissue %!in% c("NA"))

CNTN4_pubdata4=subset(CNTN4_pubdata3, seurat_clusters %!in% c(3,5,7,9,10,14,15))

sample_gt30_4=names(table(CNTN4_pubdata4$sampleID)[table(CNTN4_pubdata4$sampleID)>=10])   ### 117 samples and 11709 cells
CNTN4_pubdata4=subset(CNTN4_pubdata4, sampleID %in% c(sample_gt30_4))

DimPlot(CNTN4_pubdata4, group.by = "seurat_clusters", reduction = "umap",label = T)
CNTN4_pubdata4 <- RunUMAP(object = CNTN4_pubdata4, assay = "RNA", reduction = "harmony", metric = "cosine",
                 dims = 1:20,
                 n.neighbors = 120,
                 min.dist = 0.75, spread=1)
scPioneer::FeaturePlot2(CNTN4_pubdata4, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
table(CNTN4_pubdata4$meta_tissue)

CNTN4_pubdata3_sub <- RunUMAP(object = CNTN4_pubdata3_sub, assay = "RNA", reduction = "harmony", metric = "cosine",
                          dims = 1:20,
                          n.neighbors = 120,
                          min.dist = 0.75, spread=1)
CNTN4_pubdata3_sub <- FindClusters(object = CNTN4_pubdata3_sub, algorithm = 2,
                      resolution = c(0.5),
                      verbose = F,
                      n.start = 10)

OutputDir="./AnaRes20251019/"
table(CNTN4_pubdata3_sub$meta_tissue)
DimPlot(CNTN4_pubdata3_sub, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(CNTN4_pubdata3_sub, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)

sample_gt30_x=names(table(CNTN4_pubdata3_sub$sampleID)[table(CNTN4_pubdata3_sub$sampleID)>=10])   ### 117 samples and 11709 cells
CNTN4_pubdata3_sub=subset(CNTN4_pubdata3_sub, sampleID %in% c(sample_gt30_x))

totalX=table(total_stat_data$meta_tissue)
totalY=table(CNTN4_pubdata3_sub$meta_tissue)
totalXY=intersect(names(totalX[totalX>500]),names(totalY[totalY>50]))
totalXY=totalXY[totalXY %!in% c("testis","skin")]
df_tmp=table(CNTN4_pubdata3_sub$meta_tissue)[totalXY]/totalX[totalXY]

CNTN4_pubdata3_sub0=subset(CNTN4_pubdata3_sub, meta_tissue %in%  c(totalXY)) ### 93.5%
CNTN4_pubdata3_sub0$meta_tissue=factor(CNTN4_pubdata3_sub0$meta_tissue, levels = tmp_order)
tmp_order=c("pancreas", "Esophagus","Duodenum",
             "Mid-jejunum", "Proximal-jejunum",  "Ileum", "Proximal-small-intestine", "small intestine","large intestine","Colon",
            "Lung-airway", "Lung-airway-trachea", "Lung-tracheal-epi", "AdultKidney", "AdultThyroid", "ovary")
DimPlot(CNTN4_pubdata3_sub0,
        group.by = "meta_tissue", reduction = "umap",label = T)
p=scPioneer::DimPlot_idx(CNTN4_pubdata3_sub0, group.by = "meta_tissue", repel = T, reduction = "umap",
                       cols=berryFunctions::addAlpha(hue_pal()(16),0.3))+NoAxes()
runFig(paste0(OutputDir,"/", paste0("Fig3.pubdata_CNTN4_tissue")),p,6.5,4.5)

p=scPioneer::DimPlot_idx(CNTN4_pubdata3_sub0, group.by = "meta_tissue", ncol=4, label.idx.size = 0,
                         split.by = "meta_tissue",  reduction = "umap",
                         cols=berryFunctions::addAlpha(hue_pal()(16),0.3))
runFig(paste0(OutputDir,"/", paste0("Fig3.pubdata_CNTN4_tissue_split")),p,10,8)

p=scPioneer::FeaturePlot2(CNTN4_pubdata3_sub0,
                        features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"),
                        ncol = 4, order = T)
runFig(paste0(OutputDir,"/", paste0("Fig3.pubdata_CNTN4_marker_show")),p,10,4.8)

CNTN4_pubdata3_sub00=reDRClusterTest(CNTN4_pubdata3_sub00, sp="hum", ndim=20, nClu=100, mdist=0.25,
                              formalFeature=NULL, mypca="harmony", harmony=10, prepro=TRUE,
                              mybatch="sampleID", umap=TRUE, vstN=2000, regVar=NULL, rmvHVGribo=FALSE)  ###c("stats_mitoPercent")
CNTN4_pubdata3_sub00 <- RunUMAP(object = CNTN4_pubdata3_sub00, assay = "RNA", reduction = "harmony", metric = "cosine",
                              dims = 1:20,
                              n.neighbors = 120,
                              min.dist = 0.75, spread=1)
CNTN4_pubdata3_sub00 <- FindClusters(object = CNTN4_pubdata3_sub00, algorithm = 2,
                                   resolution = c(0.5),
                                   verbose = F,
                                   n.start = 10)

scPioneer::FeaturePlot2(CNTN4_pubdata3_sub00, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
DimPlot(CNTN4_pubdata3_sub00, sizes.highlight = 0.1, cells.highlight = colnames(subset(CNTN4_pubdata3_sub, meta_tissue %in% c("Ileum"))))
CNTN4_pubdata3_sub00$mark=0
tmpvec=CNTN4_pubdata3_sub00$mark
tmpvec[names(tmpvec) %in% colnames(subset(CNTN4_pubdata3_sub, meta_tissue %in% c("Ileum")))]=1
CNTN4_pubdata3_sub00$mark=as.vector(tmpvec)
table(CNTN4_pubdata3_sub00$seurat_clusters,CNTN4_pubdata3_sub00$mark) %>% rowPercents()

tmp_color=berryFunctions::addAlpha(hue_pal()(16),0.5)
names(tmp_color)=tmp_order
df_tmp_X=data.frame(tissue=names(df_tmp), frac=as.vector(df_tmp))
df_tmp_X=df_tmp_X[order(df_tmp_X$frac, decreasing = T),]
df_tmp_X$tissue=factor(df_tmp_X$tissue, levels = unique(df_tmp_X$tissue))
p=ggplot(df_tmp_X)+geom_bar(aes(tissue, frac), stat="identity", fill=tmp_color, color="white")+theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
runFig(paste0(OutputDir,"/", paste0("Fig3.pubdata_CNTN4_frac")),p,6,3.5)

CNTN4_pubdata3_sub$meta_tissue=factor(CNTN4_pubdata3_sub$meta_tissue, levels = unique(df_tmp_X$tissue))
p=VlnPlot(CNTN4_pubdata3_sub, group.by = "meta_tissue", features = "CNTN4", ncol = 1, pt.size = 0)+NoLegend()

library(CytoTRACE)
xdata_results <- CytoTRACE(as.matrix(CNTN4_pubdata3_sub@assays$RNA$counts),
                           ncores = 4, subsamplesize = 1000)
ptV=CNTN4_pubdata3_sub$meta_tissue %>% as.character()
names(ptV)=names(CNTN4_pubdata3_sub$meta_tissue)
### NOTE DO NOT load scPioneer!!! ->boxplot will be replaced by this package!!!
plotCytoTRACE(xdata_results, phenotype = ptV,  outputDir = "./PubData/",
              emb = CNTN4_pubdata3_sub@reductions$umap@cell.embeddings %>% as.data.frame()) #, gene = "Kit")

saveRDS(CNTN4_pubdata, "./PubData/CNTN4_pubdata.rds")
saveRDS(CNTN4_pubdata2, "./PubData/CNTN4_pubdata2.rds")
saveRDS(CNTN4_pubdata3, "./PubData/CNTN4_pubdata3.rds")
saveRDS(CNTN4_pubdata3_sub, "./PubData/CNTN4_pubdata3_sub.rds")
saveRDS(CNTN4_pubdata4, "./PubData/CNTN4_pubdata4.rds")

#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
pub_data=readRDS("./data/PubData/Immune_final_pred.rds")
pub_data$cb=rownames(pub_data@meta.data)
dim(pub_data)

data_genes=read.csv("./PubData/Immune_feature.csv")
mygenes=intersect(pcgenes$V1, data_genes$features)

pub_data_tmp=subset(pub_data, label %in%
                     c("Mono01_CD14.c","Mono02_FCGR3A","Mono03_CD14.int",
                       "Mph01_FABP4","Mph02_RGS1","Mph03_FOLR2","Mph04_MT1G","Mph05_CCL20"))
pub_data_STARD13=subset(pub_data, STARD13>0)
hhh2=c(colnames(pub_data_STARD13), colnames(pub_data_tmp)) %>% unique()

Mph_pubdata=subset(pub_data, cb %in% hhh2)
Mph_pubdata=subset(Mph_pubdata, features = mygenes)
Mph_pubdata=UpdateSeuratObject(Mph_pubdata)
Mph_pubdata=subset(Mph_pubdata, nFeature_RNA>=200)

sample_gt30=names(table(Mph_pubdata$sampleID)[table(Mph_pubdata$sampleID)>=10])   ### 669 samples and 11709 cells
Mph_pubdata=subset(Mph_pubdata, sampleID %in% c(sample_gt30))

Mph_pubdata=reDRClusterTest(Mph_pubdata, sp="hum", ndim=20, nClu=100, mdist=0.25,
                            formalFeature=NULL, mypca="harmony", harmony=10, prepro=TRUE,
                            mybatch="sampleID", umap=TRUE, vstN=2000, regVar=NULL, rmvHVGribo=FALSE)  ###c("stats_mitoPercent")

DimPlot(Mph_pubdata, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(Mph_pubdata, pt.size = 1,
                        features = c("PTPRC","STARD13","FOLR2","RGS1","FCGR3A","FABP4",
                                     "CD3D","CD79A","CPA3","CSF3R","NKG7","IL4I1"), ncol = 4, order = T)
scPioneer::FeaturePlot2(Mph_pubdata, pt.size = 1,
                        features = c("CD68","CD14","CD163","FOLR2","TIMD4","LYVE1","CCR2","MRC1","FCN1", "VCAN","CST3","CPVL"), ncol = 4, order = T)
scPioneer::DotPlot2(Mph_pubdata,  features =
                      unique(c("PTPRC","STARD13","FOLR2","RGS1","FCGR3A","FABP4",
                          "CD3D","CD79A","CPA3","CSF3R","NKG7","IL4I1",
                        "CD68","CD14","CD163","FOLR2","TIMD4","LYVE1","CCR2","MRC1","FCN1", "VCAN","CST3","CPVL")))
scPioneer::DimPlot_idx(Mph_pubdata,  group.by = "seurat_clusters", reduction = "umap")

Idents(Mph_pubdata)="seurat_clusters"
Mph_pubdata_degs=FindAllMarkers(Mph_pubdata, group.by = "seurat_clusters", only.pos = T, max.cells.per.ident = 2000, test.use = "t")

DotPlot(Mph_pubdata, features = (Mph_pubdata_degs %>% group_by(cluster) %>% top_n(6, wt =avg_log2FC))$gene %>% unique(),
        group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")

scPioneer::DimPlot_idx(subset(Mph_pubdata, label %in% c("Mono01_CD14.c","Mono02_FCGR3A","Mono03_CD14.int",
                                                        "Mph01_FABP4","Mph02_RGS1","Mph03_FOLR2","Mph04_MT1G","Mph05_CCL20")),
                       group.by = "label", reduction = "umap")

xxx=subset(Mph_pubdata, label %in% c("Mph02_RGS1","Mph03_FOLR2"))
xxx_gt30=names(table(xxx$sampleID)[table(xxx$sampleID)>=10])   ### 344 samples and 51797 cells
xxx=subset(xxx, sampleID %in% c(xxx_gt30))

Mph_pubdata_subTotal=subset(subset(Mph_pubdata, label %in% c("Mono01_CD14.c","Mono02_FCGR3A","Mono03_CD14.int",
                                                             "Mph01_FABP4","Mph02_RGS1","Mph03_FOLR2","Mph04_MT1G","Mph05_CCL20")),
                            sampleID %in% xxx_gt30)
totalX=table(Mph_pubdata_subTotal$meta_tissue)
df_tmp=table(xxx$meta_tissue)[names(totalX)]/totalX
df_tmp_X=data.frame(tissue=names(df_tmp), frac=as.vector(df_tmp))
df_tmp_X=df_tmp_X[order(df_tmp_X$frac, decreasing = T),]
df_tmp_X$tissue=factor(df_tmp_X$tissue, levels = unique(df_tmp_X$tissue))
ggplot(df_tmp_X)+geom_bar(aes(tissue, frac), stat="identity", color="grey")+theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

xxx$meta_tissue=factor(xxx$meta_tissue, levels = unique(df_tmp_X$tissue))
VlnPlot(xxx, group.by = "meta_tissue", features = c("RGS1","FOLR2","STARD13"),
        ncol = 1, pt.size = 0, stack = T, flip = T)+NoLegend()

FeaturePlot(subset(xxx, meta_tissue %in% c(head(df_tmp_X, n=3)$tissue %>% as.vector())),
            features = c("RGS1","FOLR2","STARD13"), order = T, ncol = 3)

plot_list=list()
for(xi in df_tmp_X[22:23,]$tissue %>% as.vector()){
  ppp=FeatureScatter(subset(xxx, meta_tissue %in% c(xi)),
                                 feature1 = "STARD13", feature2 = "RGS1",
                                 group.by = "label")
  ppp_rate=paste0(paste(c("Neg","sPos","dPos"),table(rowSums(ppp$data[,1:2]>0))),collapse = "; ")
  qqq=FeatureScatter(subset(xxx, meta_tissue %in% c(xi)),
                                 feature1 = "STARD13", feature2 = "FOLR2",
                                 group.by = "label")
  qqq_rate=paste0(paste(c("Neg","sPos","dPos"),table(rowSums(qqq$data[,1:2]>0))),collapse = "; ")
  plot_list[[xi]]=
    ppp+ggtitle(ppp_rate)+theme(plot.title = element_text(size=5))+
    qqq+ggtitle(qqq_rate)+theme(plot.title = element_text(size=5))
}
plot_grid(plot_list[[1]],plot_list[[2]],#plot_list[[3]],plot_list[[4]],plot_list[[5]],
          ncol = 1, labels = names(plot_list), label_size = 10)

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]], plot_list[[6]],
          ncol = 2, labels = names(plot_list), label_size = 10)

#######################################
#######################################
#######################################
load("./PubData/test.image")

OKOK=read.csv2("./PubData/Immune_barcode.csv")

#######################################
#######################################
#######################################
step1_QC=function(tmpSdir, tmpSname, tmpSgroup, QCdir="./", Species="mm",
                  pMT=25, nF=200, ...){

  if (Species == "mm") {    ### CAUTION: maybe not complete lists in mouse
    MitochondrialPattern <- "^mt\\." #"^mt-"
    RiboPattern <- "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa"
    IGPattern <- "^Ig[klh].+"
    HBPattern <- "^Hb[abq].+" #"^Hb[abq][1-].+"
  }else{
    MitochondrialPattern <- "^MT-"
    RiboPattern <- "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
    IGPattern <- "^IG[KLH].+"
    HBPattern <- "^HB[ABDEGMQZ].?$"
  }

  tmpSdir="./"
  file=list.files(tmpSdir, recursive = F, pattern = "*\\.gct$")
  new_counts <- read.table(file=paste0(tmpSdir,"/",file[1]), sep = "\t", header = TRUE)
  geneidx=new_counts$Description
  geneidx_rep=table(geneidx)[table(geneidx)>1] %>% names()
  new_counts=new_counts[new_counts$Description %!in% geneidx_rep,]
  cellidx=colnames(new_counts)[3:length(colnames(new_counts))]
  new_counts=new_counts[,cellidx]
  colnames(new_counts)=gsub("\\.","_",colnames(new_counts))
  #  geneidx[match(tmpGdf[i,]$id1,geneidx)]=tmpGdf[i,]$id2
  #}
  rownames(new_counts)=geneidx[geneidx %!in% geneidx_rep]
  obj_xxx <- CreateSeuratObject(counts = new_counts, min.cells = 3, project = "GTExV8")
  obj_xxx$cellid=colnames(obj_xxx)
  seurat_obj=obj_xxx
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = MitochondrialPattern)
  seurat_obj<- PercentageFeatureSet(seurat_obj, pattern = RiboPattern, col.name = "percent.ribo")
  seurat_obj<- PercentageFeatureSet(seurat_obj, pattern = IGPattern, col.name = "percent.igklh")
  seurat_obj<- PercentageFeatureSet(seurat_obj, pattern = HBPattern, col.name = "percent.hbb")
  meta_add=read.table("./GTEx_Analysis_v8_RNAseq_samples.txt",header = T)
  meta_add$cellid=gsub("-","_",meta_add$sample_id)
  oldmeta=seurat_obj@meta.data
  oldmeta=left_join(oldmeta, meta_add, by="cellid")
  rownames(oldmeta)=oldmeta$cellid
  seurat_obj=AddMetaData(seurat_obj, oldmeta)

  xfun=function(mData){
    UMICellCount=data.frame(as.vector(table(mData$orig.ident)), row.names = names(table(mData$orig.ident)))
    UMICountMedian <- as.data.frame(tapply(mData$nCount_RNA, INDEX=mData$orig.ident, FUN=median))
    FeatureMedian <- as.data.frame(tapply(mData$nFeature_RNA, INDEX=mData$orig.ident, FUN=median))
    MitochondrialPercentageMedian <- as.data.frame(tapply(mData$percent.mt, INDEX=mData$orig.ident, FUN=median)) %>% round(digits = 2)
    RiboPercentageMedian <- as.data.frame(tapply(mData$percent.ribo, INDEX=mData$orig.ident, FUN=median)) %>% round(digits = 2)
    IgPercentageMedian <- as.data.frame(tapply(mData$percent.igklh, INDEX=mData$orig.ident, FUN=median)) %>% round(digits = 2)
    HbbPercentageMedian <- as.data.frame(tapply(mData$percent.hbb, INDEX=mData$orig.ident, FUN=median)) %>% round(digits = 2)
    print(paste0(UMICellCount," | ",UMICountMedian," | ",FeatureMedian," | ",
                 MitochondrialPercentageMedian," | ",RiboPercentageMedian," | ",IgPercentageMedian," | ",HbbPercentageMedian))
  }
  xfun(seurat_obj)
  seurat_obj=seuQuickRun(seurat_obj, batchsize = 1)

}

DimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "umap",label = T)
DimPlot(seurat_obj, group.by = "tissue_id", reduction = "umap",label = T)

MP_resident_markers=c("FOLR2","TIMD4","LYVE1","CD74","CCR2","MRC1", "CD163", "SIGLEC1","IGF1","PDGFC","CD68","CD14")
seu_o=subset(Mobj_keep, MajorClass02 %in% c("Epithelial Cells"))
seu_o$subtypes015=droplevels(seu_o$subtypes015)
xxx_markers=FindMarkers(seu_o, ident.1 = "CNTN4_NEL_Cells", group.by = "subtypes015", logfc.threshold = 0.5, test.use = "t", min.pct = 0.35, only.pos = T)
seu_e=subset(Mobj_keep, MajorClass02 %in% c("Immune Cells"))
seu_e=subset(seu_e, subtypes015 %in% c("DC01_CD1C","Mono03_CD14.int","Mph01_FABP4","Mph02_RGS1",
                                            "DC02_CD207","Mono01_CD14.c","Mono02_FCGR3A",
                                            "Mph03_FOLR2","Mph04_MT1G","Mph05_CCL20",
                                            "Neutro01_NAMPT","Neutro02_MMP9","Neutro03_CXCL8"))
seu_e$subtypes015=droplevels(seu_e$subtypes015)
yyy_markers=FindMarkers(seu_e, ident.1 = c("Mph02_RGS1","Mph03_FOLR2"), group.by = "subtypes015", logfc.threshold = 0.5, test.use = "t", min.pct = 0.35, only.pos = T)

seurat_obj=AddModuleScore(seurat_obj,
                          features = list(rownames(yyy_markers)),
                          name = "ResidMph_score",
                          ctrl = 100)
seurat_obj=AddModuleScore(seurat_obj,
                          features = list(rownames(xxx_markers)),
                          name = "CNTN4_score",
                          ctrl = 100)
p0=VlnPlot(seurat_obj, features = c("ResidMph_score1","CNTN4_score1"),
        group.by = "tissue_id", ncol = 1, pt.size = 0)
vln_modif(p0)

ggplot(seurat_obj@meta.data)+
  geom_point(aes(ResidMph_score1,CNTN4_score1,color=tissue_id))+ggtitle(paste0("Cor:","0.633"))

seurat_obj=AddModuleScore(seurat_obj,
                          features = list(setdiff(rownames(yyy_markers),
                                                  intersect(rownames(yyy_markers), rownames(xxx_markers)))),
                          name = "ResidMph_NOverlap_score",
                          ctrl = 100)
seurat_obj=AddModuleScore(seurat_obj,
                          features = list(setdiff(rownames(xxx_markers),
                                                   intersect(rownames(yyy_markers), rownames(xxx_markers)))),
                          name = "CNTN4_NOverlap_score",
                          ctrl = 100)
seurat_obj=readRDS(paste0("./", "GTEx_Analysis_v8", ".01.RawSeu.rds"))
p=ggplot(seurat_obj@meta.data)+
  geom_point(aes(ResidMph_NOverlap_score1,CNTN4_NOverlap_score1,color=tissue_id))+ #+ggtitle(paste0("Cor:","0.633"))
  geom_smooth(aes(ResidMph_NOverlap_score1,CNTN4_NOverlap_score1), method=lm,se = T)+
  theme_cowplot()
runFig(paste0(OutputDir,"/", paste0("Fig4.CNTN4_residMph_cor","_GTEx")), p, 14, 4.5)
plot_list=list()
for(x in unique(seurat_obj@meta.data$tissue_id)){
  tmp_df=seurat_obj@meta.data[seurat_obj@meta.data$tissue_id %in% x,]
  pxx=cor.test(tmp_df[,"ResidMph_NOverlap_score1"],tmp_df[,"CNTN4_NOverlap_score1"])
  p=ggplot(tmp_df)+
    geom_point(aes(ResidMph_NOverlap_score1,CNTN4_NOverlap_score1),color="darkgrey")+
    ggtitle(paste0(x,"\nCor:",round(as.vector(pxx$estimate),digits = 3),"; Pval:",round(pxx$p.value, digits = 4)))+
    geom_smooth(aes(ResidMph_NOverlap_score1,CNTN4_NOverlap_score1), method=lm,se = T)+
    theme_cowplot()
  plot_list[[x]]=p
}
runFig(paste0(OutputDir,"/", paste0("Fig4.LARGE_CNTN4_residMph_cor","_GTEx_splitbyTissue")),
       plot_grid(plotlist = plot_list, ncol = 5), 14, 50)
cor.test(seurat_obj@meta.data[,"ResidMph_NOverlap_score1"],seurat_obj@meta.data[,"CNTN4_NOverlap_score1"])
### cor.test: 0.438 (after removing 163 overlapped genes)

CNTN4_data_ori=AddModuleScore(CNTN4_data_ori,
                              features = list(rownames(xxx_markers)),
                              name = "CNTN4_1score",
                              ctrl = 100)
CNTN4_data_ori=AddModuleScore(CNTN4_data_ori,
                              features = list(c("CNTN4","GPC6","FBXL7","YAP1")),
                              name = "CNTN4_2score",
                              ctrl = 100)
ggplot(CNTN4_data_ori@meta.data)+
  geom_point(aes(CNTN4_1score1,neuron1,color=types))+ggtitle(paste0("Cor:","0.206"))+
ggplot(CNTN4_data_ori@meta.data)+
  geom_point(aes(CNTN4_2score1,neuron1,color=types))+ggtitle(paste0("Cor:","0.264"))
cor.test(CNTN4_data_ori@meta.data[,"CNTN4_2score1"],CNTN4_data_ori@meta.data[,"neuron1"])

vln_modif=function(obj,...){
  for (i in 1:length(obj$patches$plots)){
    obj$patches$plots[[i]]$labels$x=NULL
    obj$patches$plots[[i]]=obj$patches$plots[[i]]+
      geom_boxplot(width=0.35)+
      theme(axis.text.x = element_blank ())
  }
  obj$labels$x=NULL
  return(obj+geom_boxplot(width=0.35))
}
seuQuickRun=function(seu,batchsize=1){
  seu = seu[,sample(colnames(seu),
                    round(batchsize*dim(seu)[2]))] %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose=FALSE) %>%
    RunUMAP(reduction="pca", dims=1:20) %>%
    FindNeighbors(reduction="pca", dims=1:20, verbose = FALSE) %>%
    FindClusters(resolution = c(0.5), verbose = FALSE)
  Idents(seu)="seurat_clusters"
  return(seu)
}

devdata01 <- Read10X(data.dir = "./devData/")
seu_devdata01 = CreateSeuratObject(counts = devdata01)
devdata01_addmeta <- read.table("./devData/barcodes_meta.tsv", header = T, sep = "\t")
rownames(devdata01_addmeta)=devdata01_addmeta$cell_id
seu_devdata01 <- AddMetaData(seu_devdata01, devdata01_addmeta)

seu_devdata01=seuQuickRun(seu_devdata01,batchsize = 1)

#######################################
#######################################
#######################################
xxx=readRDS("DATA_DIR")
xxx$Celltype=xxx$MajorCellTypes0701
xxx_anno=xxx@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","Sample","cellid","Tissue","Species","Group","Celltype")]
xxx_umap=xxx@reductions$umap@cell.embeddings
tmpV=xxx_umap[,2]
xxx_umap[,2]=xxx_umap[,1]
xxx_umap[,1]=tmpV
write.csv(xxx_umap, file="DATA_DIR")
write.csv(xxx_anno, file="DATA_DIR")

Blike03_CtlAct=readRDS("DATA_DIR")
indata=Blike03_CtlAct

tmpdf=Blike03_CtlAct@meta.data[,c("cellid","Bsubtypes")]
tmpdf$Bsubtypes=as.character(tmpdf$Bsubtypes)
tmpdf$Bsubtypes[tmpdf$cellid %in% colnames(subset(memBdata, seurat_clusters %in% c(4,5,6,15)))]="Switched memory B cells"   ##c(7,8,10,13)
#tmpdf$Btypes03[tmpdf$cellid %in% colnames(subset(memBdata, seurat_clusters %in% c(16)))]="IGHG_smB"  ## only 24 cells
Blike03_CtlAct$Bsubtypes=tmpdf$Bsubtypes
Blike03_CtlAct$Bsubtypes=factor(Blike03_CtlAct$Bsubtypes, levels = c("Naive B cells","Memory B cells","Switched memory B cells","Plasma"))

Idents(Blike03_CtlAct)="Bsubtypes"
DimPlot(Blike03_CtlAct, reduction = "umap", label = T, repel = T, raster = T) + NoLegend()

###---

tmpdf=Blike03_RemRla@meta.data[,c("cellid","Bsubtypes")]
tmpdf$Bsubtypes=as.character(tmpdf$Bsubtypes)
tmpdf$Bsubtypes[tmpdf$cellid %in% colnames(subset(memBdata, seurat_clusters %in% c(4,5,6,15)))]="Switched memory B cells"   ##c(7,8,10,13)
#tmpdf$Btypes03[tmpdf$cellid %in% colnames(subset(memBdata, seurat_clusters %in% c(16)))]="IGHG_smB"  ## only 24 cells
Blike03_RemRla$Bsubtypes=tmpdf$Bsubtypes
Blike03_RemRla$Bsubtypes=factor(Blike03_RemRla$Bsubtypes, levels = c("Naive B cells","Memory B cells","Switched memory B cells","Plasma"))

Idents(Blike03_RemRla)="Bsubtypes"
DimPlot(Blike03_RemRla, reduction = "umap", label = T, repel = T, raster = T) + NoLegend()

Seurat::DimPlot(indata, group.by = "Bsubtypes")

################################################################################
subNeu=readRDS("DATA_DIR")

xxx=Seurat::VlnPlot(subset(subNeu, subset=TNFSF13>0), pt.size = 0, slot = "data",
            assay = "RNA", features = c("TNFSF13"), group.by = "Group", split.by = "Group")+
  NoLegend()+xlab("")+ylab("")   #geom_boxplot()+
xxx$layers[[1]]$stat_params$trim=FALSE
xxx=xxx+stat_summary(fun.y = median, geom = "point", shape=95, size = 25, color = "darkred", position = position_dodge(width = 0.9)) #geom_boxplot()
xxx=xxx+NoLegend()+xlab("")+scale_fill_manual(values = c(Control="lightblue", Active="orange"))
runFig(paste0(FinalFig, "FigureXX", ".CtlAct.APRIL.vln"), xxx, 2.5, 4.25)
t.test(xxx$data[xxx$data$ident=="Control",]$TNFSF13, xxx$data[xxx$data$ident=="Active",]$TNFSF13)

Seurat::VlnPlot(subNeu, features = "TNFSF13", group.by="NEUsubtypes")
Seurat::DotPlot(subNeu, features = "TNFSF13", group.by="Group")

#######################################################
indata=readRDS("DATA_DIR")
write.csv(table(indata$Sample,indata$BCR_ClonotypeCount), "DATA_DIR")
indata=readRDS("DATA_DIR")
write.csv(table(indata$Sample,indata$BCR_ClonotypeCount), "DATA_DIR")

indata$BCR_class=ifelse(indata$BCR_ClonotypeCount == "0",
                        yes = "None",
                        no = ifelse(indata$BCR_ClonotypeCount == "1",
                                    yes = "Unique",
                                    no = ifelse(indata$BCR_ClonotypeCount %in% c("2","3","4"),
                                                yes = as.character("≥2"),
                                                no = ifelse(indata$BCR_ClonotypeCount %in% c("5","6","7","8","9"),
                                                            yes = as.character("≥5"),
                                                            no = as.character("≥10")))))
indata=subset(indata, BCR_class %!in% c("None"))

mycolor=c(hue_pal()(4)[3], hue_pal()(12)[4], hue_pal()(12)[2], hue_pal()(4)[1])
cond=c("Control","Active")  ###,"REM","RLA")
indata$Group=factor(indata$Group, levels =cond)
tmpdf=RcmdrMisc::colPercents(table(indata[["BCR_class"]][,1], indata$Group))
tmpdf=as.data.frame(tmpdf[1:(dim(tmpdf)[1]-2),])
tmpdf$Classtype=rownames(tmpdf)
tmpdf=tmpdf[,c(cond, "Classtype")]
tmpdf=tmpdf[rev(c("Unique","≥2","≥5","≥10")),]

tmpdf=tidyr::pivot_longer(tmpdf, cols = Control:Active, names_to = "Condition", values_to = "Percentage")

donut_plot=list()
for(i in c(cond)){   ### "MS","HC"
  subdf=as.data.frame(subset(tmpdf, Condition==i))

  subdf$Classtype=factor(subdf$Classtype, levels = as.character(c("Unique","≥2","≥5","≥10")))

  subdf$ymax = cumsum(subdf$Percentage)
  subdf$ymin = c(0, head(subdf$ymax, n=-1))
  subdf$labelPosition <- (subdf$ymax + subdf$ymin) / 2
  subdf$label <- paste0(subdf$Percentage,"%")
  donut_plot[[i]]=
    ggplot(subdf, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Classtype), shape=18) +
    geom_rect(color="black") +
    geom_label( x=3.5, aes(y=labelPosition, label=label), label.size = 0, show.legend = F) + myfont +
    scale_fill_manual(values = mycolor) +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(1.5, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() +
    theme(legend.position = "top", legend.direction = "horizontal") + labs (fill='NO. of clones') +
    guides(fill = guide_legend(override.aes=list(shape = 20)))#+
}

#xxx=plot_grid(plotlist = donut_plot)
xxx=plot_grid(plotlist = donut_plot)
runFig(paste0(cdir, "Figure1h", ".CtlAct.donut_Blike.label"), xxx, 11, 6)

#xxx=plot_grid(plotlist = donut_plot)
#xxx=plot_grid(plotlist = donut_plot)

xxx=table(indata$Sample,indata$BCR_ClonotypeCount)
xxx=xxx[,c(-1,-2)]
xi=xj=c()
for(x in rownames(xxx)[1:16]){
  xline=xxx[x,]
  xline=as.vector(xline)
  if(length(xline)>0){
    iSimpson=1/sum((xline/sum(xline))^2)
  }else{
    iSimpson=0
  }
  xi=c(xi,x)
  xj=c(xj,iSimpson)
}
ctlact=data.frame(sample=xi, iSimpson=xj)

#xxx@reductions$umap_bak2=xxx@reductions$umap
xxx=RunUMAP(xxx, dims = 1:16, n.neighbors = 350, spread = 1.75, #repulsion.strength = 1.5,
            min.dist = 0.05, seed.use = 1, n.epochs = 300)
DimPlot(xxx, group.by = "major_type_bak", label = T)
#xxx$major_type_bak=xxx$major_type
#xxx=modifyAno(xxx, anoCol = "major_type", oldID = "Pericyte", newID = "FIB")
