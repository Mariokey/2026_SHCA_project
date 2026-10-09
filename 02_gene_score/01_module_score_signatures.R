# ============================================================================
# Signature-based module scores with AddModuleScore (monocyte / macrophage states)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Gene score analysis"
# Source     : extracted from Analysis240525.r (L6019-6248, L8151-8277)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
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

dataSubHSC=subset(HSPCdata, (types01=='HSC'))
STEMsigG_old=STEMgenes
CMPsigG_old=unique(c("PRSS57", "CDK6",   "HSPD1",  "LGALS1", "CEBPA",
                     CMP_degs[1:50,]$gene))
GMPsigGenes_old=unique(c("AZU1",  "MPO",   "CTSG",  "ELANE", "LYST", "CAMP", "LCN2", "RETN", "LTF",  "MNDA",
                         "LYZ", "CST3", "RETN", "FCER1G", "TYROBP", "S100A8", "ANXA2", "FTL", "GRN", "S100A9", "RAB31", "PRTN3", "STXBP2", "MNDA", "CTSZ", "CTSG", "MPO", "AZU1", "PLAC8", "CFD", "TUBB4B", "S100A11", "PRDX4", "ASAH1",
                         "LYZ", "FCER1G", "S100A8", "RNASE6", "GRN", "S100A9", "IRF8", "CTSG", "MPO", "AZU1", "PLAC8", "ISG15", "HMGB2"))
LYMsig1_old=unique(c(formalFeature[["CLP"]], formalFeature[["BPro"]],formalFeature[["BPre"]], "MS4A1", "CD79A",
                     CLP_degs[1:50,]$gene))

out0203="DATA_DIR"
dir.create(out0203)
pptTEMP = paste0("DATA_DIR")
my_ppt <- read_pptx(path = pptTEMP)
outPPT = paste0(out0203, "/tjLHD_V20230203_Heatmap_GeneModule_03.CtlAct.pptx")
for (p in seq(0.7,0.95,0.05)){
  runHeat(dataSubHSC, qt=c(p,p+0.05), targetCT=c("HSC", "CMP", "GMP", "CLP"), type=TRUE)
}
print(my_ppt, target = outPPT)
rm(my_ppt)

### Immu & MP refine 1113
subdataList=readRDS("./subdataList_241011.rds")
ImmuLData=subdataList$lym
ImmuMData=subdataList$myl
ImmuLData=subset(ImmuLData,SubClass00 %!in% c("CD8T01_hsp"))
ImmuLData$SubClass00=droplevels(ImmuLData$SubClass00)
DimPlot_idx(ImmuLData, group.by = "SubClass00", cols=scPalette2(17), prefix.index = 1:17)
ImmuMData$SubClass01=ImmuMData$SubClass00
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mono03_THBS1", newID = "Mono02_CD14")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mono01_CD14", newID = "Mono01_FCGR3A")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mono03_VCAN", newID = "Mono01_FCGR3A")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mono04_LYPD2", newID = "Mono01_FCGR3A")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mono02_INHBA", newID = "Mph03_INHBA")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mph01_AGRP", newID = "Mph01_FABP4")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mph04_CD5L", newID = "Mph02_STARD13")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mph05_SIGLEC8", newID = "Mph01_FABP4")
ImmuMData=modifyAno(ImmuMData, anoCol = "SubClass01", oldID = "Mph03_HLA", newID = "DC01_CD1C")
DimPlot_idx(ImmuMData, group.by = "SubClass01", cols=scPalette2(9), prefix.index = 18:26)

plot11=roveRun(ImmuLData, Group="Location_short", CT="SubClass00", outdir=paste0(dir_img,"L."), flip = FALSE, mywidth=4.5, nCT=17)
plot12=roveRun(ImmuMData, Group="Location_short", CT="SubClass01", outdir=paste0(dir_img,"M."), flip = FALSE, mywidth=4.5, nCT=9)

xxx1=read.table(paste0(dir_img,"L.rove_main_age.txt"))
xxx2=read.table(paste0(dir_img,"M.rove_main_age.txt"))
write.table(rbind(xxx1,xxx2), sep = '\t', quote = F, file = paste0(dir_img,'LM.rove_main_age.txt'))
ImmuLData$Location_short=factor(ImmuLData$Location_short, levels = loca_order)
ImmuMData$Location_short=factor(ImmuMData$Location_short, levels = loca_order)
xxx1=table(ImmuLData@meta.data[,'SubClass00'],ImmuLData$Location_short)
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
FeDM_score=c("FOLR2","TIMD4","RGS1","LYVE1","NINJ1","GAS6","MRC1","CD163","IGFB4","MERTK","ABCA1","SIGLEC1","CSF1R","SELENOP")

xdata=mp_trajec
xdata[["RNA"]]=as(xdata[["RNA"]], Class = "Assay")
Idents(xdata)="seurat_clusters"
xdata=angryAno(xdata, outdir = dir_img)
xdata[["RNA"]]=as(xdata[["RNA"]], Class = "Assay5")

DotPlot(xdata, features = c("FCN1","STARD13","MCEMP1","FABP4","INHBA","APOC1",MoDM_score), group.by = "Annotation_angrycell")+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")

Idents(xdata)="seurat_clusters"
xxx=FindAllMarkers(xdata, group.by = "seurat_clusters", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)

xdata=subset(xdata, seurat_clusters %!in% c(6,7,8))
xdata$seurat_clusters=droplevels(xdata$seurat_clusters)
DotPlot(xdata, features = c("FCN1","FCGR3A","CD14","STARD13","FOLR2","RGS1","FABP4","MT1G","CCL20"), group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")

myAnnList=list()
myAnnList[["MPsub"]]=c("c1_CD14+.cMo", #0
                       "c6_FABP4+.Mph",#1
                       "c4_RGS1+.Mph", #2
                       "c5_FOLR2+.Mph", #3
                       "c6_FABP4+.Mph", #4
                       "c2_CD14+.ncMo", #5
                       "c7_MT1G+.Mph", #9
                       "c8_CCL20+.Mph", #10
                       "c3_FCGR3A+.inMo") #11
Idents(xdata)="seurat_clusters"
xdata=BatchAnno_clu2ct(dataName="MPsub", xdata, annolist=myAnnList, anno_order=NULL, colnameID=paste0("MP","subtype00"))

xdata=AddModuleScore(xdata, features = list(Mono_score),
                     ctrl = 100, name = "Mono_score.0")
xdata=AddModuleScore(xdata, features = list(MoDM_score),
                     ctrl = 100, name = "MoDM_score.0")
xdata=AddModuleScore(xdata, features = list(FeDM_score),
                     ctrl = 100, name = "FeDM_score.0")

myAnnList=list()
myAnnList[["MYLsub"]]=c("Mono01_CD14.c", #0
                       "Mph01_FABP4",#1
                       "Mph02_RGS1", #2
                       "Mph03_FOLR2", #3
                       "Mph01_FABP4", #4
                       "Mono02_FCGR3A", #5: nc
                       "DC01_CD1C",
                       "RMV",
                       "RMV",
                       "Mph04_MT1G", #9
                       "Mph05_CCL20", #10
                       "Mono03_CD14.int") #11
Idents(mp_trajec)="seurat_clusters"
mp_trajec=BatchAnno_clu2ct(dataName="MYLsub", mp_trajec, annolist=myAnnList, anno_order=NULL, colnameID=paste0("MYL","subtype00"))

mp_trajec$MYLsubtype01=mp_trajec$MYLsubtype00
Idents(mp_trajec)="MYLsubtype01"
mp_trajec=modifyAno(mp_trajec, anoCol="MYLsubtype01",
                    cellid = colnames(subset(mp_trajec, MYL_SubClass05 %in% c("DC02_CD207"))),
                    newID = "DC02_CD207",
                    newOrder = c("DC01_CD1C","DC02_CD207",
                                 "Mono01_CD14.c","Mono02_FCGR3A","Mono03_CD14.int",
                                 "Mph01_FABP4","Mph02_RGS1","Mph03_FOLR2","Mph04_MT1G",
                                 "Mph05_CCL20","RMV"))
DotPlot(mp_trajec, features = c("CD1C","CD207","FCN1","FCGR3A","CD14","STARD13","FOLR2","RGS1","FABP4","MT1G","CCL20"), group.by = "MYLsubtype01")

xdata_results <- CytoTRACE(as.matrix(xdata@assays$RNA$counts),
                                          ncores = 4, subsamplesize = 1000)
ptV=xdata$MPsubtype00 %>% as.character()
names(ptV)=names(xdata$MPsubtype00)
