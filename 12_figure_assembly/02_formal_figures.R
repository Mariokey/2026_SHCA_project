# ============================================================================
# Curated formal figure panels
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Statistics and reproducibility"
# Source     : extracted from Analysis240525.r (L9349-10063)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
outdir="./AnaRes20250101/"
dir.create(outdir)

Marker_Lev1=c("KRT1","KRTDAP","DMKN",
              "PTPRC","CD3E","CD68",
              "ACTA2","MYH11","DCN",
              "CFD","CDH5","PECAM1",
              "TNP1","PRM1","DDX4",
              "MKI67","TOP2A","UBE2C")
scPioneer::FeaturePlot2(Mobj, features = Marker_Lev1, reduction = "umap",label = F, pt.size = 1, ncol = 3)
ggsave(filename = paste0(outdir,"/Mobj_Marker_Lev1.pdf"),width = 8, height = 12)
target_var=c("Donor","Germ layer","Germ layer second","Germ layer 3rd","System","Organ","Location_short")
#Donor
scPioneer::DimPlot_idx(Mobj, group.by = "Donor", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F, prefix.index = "", idx.sep = "")+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Donor.pdf"),width = 7, height = 6)
#Gender ???
#DevLayer1
scPioneer::DimPlot_idx(Mobj, group.by = "Germ layer", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Germlayer1.pdf"),width = 7, height = 6)
#DevLayer2
scPioneer::DimPlot_idx(Mobj, group.by = "Germ layer second", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Germlayer2.pdf"),width = 7, height = 6)
#DevLayer3
scPioneer::DimPlot_idx(Mobj, group.by = "Germ layer 3rd", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Germlayer3.pdf"),width = 7, height = 6)
#Sys
Mobj$System[is.na(Mobj$System)]="SKIN"
Mobj$System=factor(Mobj$System, levels = c("DS","RS","US","RPS","CS","SKIN"))
scPioneer::DimPlot_idx(Mobj, group.by = "System", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_System.pdf"),width = 7, height = 6)
#Organ
scPioneer::DimPlot_idx(Mobj, group.by = "Organ", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Organ.pdf"),width = 7, height = 6)
#Loc
scPioneer::DimPlot_idx(Mobj, group.by = "Location_short", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_Location_short.pdf"),width = 9, height = 6)

seu_epi_old=readRDS("./seu_epi_sample5w.withXepi_unk.v1216.rds")

subobjlist=list()
subobjlist["seu_epi_curat"]=list(c("forCombine01"))
toRMVcb=setdiff(as.vector(seu_epi_old$cellid),as.vector(seu_epi_curat$cellid))
Mobj_keep=subset(Mobj, cellid %!in% toRMVcb)

Mobj_keep$subtypes013=Mobj_keep$subtypes012
Mobj_keep$subtypes013=as.character(Mobj_keep$subtypes013)
Mobj_keep=addNewMetaColumn(Mobj_keep, subobjlist, nameV=NULL, forceAdd=TRUE, forceAddname="subtypes013")

sub107order=levels(Mobj_keep$subtypes012)
sub107order=c(sub107order[1:31],
              levels(seu_epi_curat$forCombine01),
              sub107order[62:107])
Mobj_keep$subtypes013=factor(Mobj_keep$subtypes013, levels = sub107order)
Mobj_keep$subtypes013=droplevels(Mobj_keep$subtypes013)

rename_list=list()
rename_list["Mono"]=list(c("Mono01_CD14.c","Mono02_FCGR3A","Mono03_CD14.int"))
rename_list["DC"]=list(c("DC01_CD1C","DC02_CD207"))
rename_list["Mph"]=list(c( "Mph01_FABP4", "Mph02_RGS1", "Mph03_FOLR2", "Mph04_MT1G", "Mph05_CCL20"))
rename_list["Intestinal"]=list(c("Enterocytes"))
rename_list["Skin"]=list(c("Sweet_Gland_Cells"))
rename_list["Umbrella"]=list(c("Umbrella_Cells"))
Mobj_keep$refine_ct02=Mobj_keep$refine_ct01
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct02", oldID = "Epi_undef1", newID = "CNTN4_Cells")
for(x in names(rename_list)){
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct02",
               cellid = colnames(subset(Mobj_keep, subtypes013 %in% rename_list[[x]])),
               newID = x)
}
sub31order=levels(Mobj_keep$refine_ct01)
sub31order=c(sub31order[1:6], "DC", sub31order[7:14],"Umbrella", sub31order[15:23], "CNTN4", sub31order[24:32])
Mobj_keep$refine_ct02=factor(Mobj_keep$refine_ct02, levels = sub31order)
Mobj_keep$refine_ct02=droplevels(Mobj_keep$refine_ct02)
write.csv(table(Mobj_keep$refine_ct02,Mobj_keep$subtypes013), "./yyy.csv")

scPioneer::DimPlot_idx(Mobj_keep, group.by = "MajorClass02", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_cellclass.pdf"),width = 7, height = 6)
scPioneer::DimPlot_idx(Mobj_keep, group.by = "refine_ct02", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_celltype.pdf"),width = 9, height = 6)
scPioneer::DimPlot_idx(Mobj_keep, group.by = "subtypes013", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_subtype.pdf"),width = 18, height = 6)

pheatmap(colPercents(table(Mobj_keep$subtypes013, Mobj_keep$refine_ct02))[1:107,], cluster_rows = F, cluster_cols = F)
ggsave(filename = paste0(outdir,"/Mobj_subtypeXcelltype.pdf"),width = 8, height = 21)

pheatmap(colPercents(table(Mobj_keep$subtypes013, Mobj_keep$Location_short))[1:107,], cluster_rows = F, cluster_cols = F)

###
outdir="./AnaRes20250110/"
dir.create(outdir)
###??? Umbrella cells和 papiliary epithelial cells合并叫urinary tract cell? BUT papiliary also in RS system (250106)
Mobj_keep$refine_ct03=Mobj_keep$refine_ct02
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct03", oldID = "^CNTN4", newID = "NE_like",
                    newOrder = gsub("^CNTN4","NE_like",levels(Mobj_keep$refine_ct03)))
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct03", oldID = "^CCDPC", newID = "Kidney",
                    newOrder = gsub("^CCDPC","Kidney",levels(Mobj_keep$refine_ct03)))
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct03", oldID = "^Papillary", newID = "Urinary_tract",
                    newOrder = gsub("^Papillary","Urinary_tract",levels(Mobj_keep$refine_ct03)))
Mobj_keep=modifyAno(Mobj_keep, anoCol="refine_ct03", oldID = "^Umbrella", newID = "Urinary_tract",
                    newOrder = levels(Mobj_keep$refine_ct03))
Mobj_keep$refine_ct03=droplevels(Mobj_keep$refine_ct03)
sub32order=c("B","Plasma","T","NK",
             "DC","Mono","Mph","Neutro","Mast",
             "Gastric","Intestinal","Endocrine","Pancreas","Hepatocyte","BEC","Tuft",
             "Basal","Superbasal","Skin","Tracheal","Urinary_tract","Kidney",
             "Granulosa","Mesothelial","NE_like",
             "EC","FB","SMC","SPG","SPC",
             "ST","Prolif.")
Mobj_keep$refine_ct03=factor(Mobj_keep$refine_ct03, levels = sub32order)

scPioneer::DimPlot_idx(Mobj_keep, group.by = "refine_ct03", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_celltype.pdf"),width = 9, height = 6)

Mobj_keep$subtypes014=Mobj_keep$subtypes013
Mobj_keep=modifyAno(Mobj_keep, anoCol="subtypes014", oldID = "^CNTN4_Cells", newID = "CNTN4_NEL_Cells",
                    newOrder = gsub("^CNTN4_Cells","CNTN4_NEL_Cells",levels(Mobj_keep$subtypes014)))
Mobj_keep=modifyAno(Mobj_keep, anoCol="subtypes014", oldID = "^MyoFB01_ELN", newID = "myoFB01_ELN",
                    newOrder = gsub("^MyoFB01_ELN","myoFB01_ELN",levels(Mobj_keep$subtypes014)))
Mobj_keep=modifyAno(Mobj_keep, anoCol="subtypes014", oldID = "^MESO01_WT1", newID = "FB14_WT1",
                    newOrder = gsub("^MESO01_WT1","FB14_WT1",levels(Mobj_keep$subtypes014)))
sub107order=c(levels(Mobj_keep$subtypes014)[1:87],
              "FB14_WT1",
              levels(Mobj_keep$subtypes014)[88:89],
              levels(Mobj_keep$subtypes014)[91:107])
Mobj_keep$subtypes014=factor(Mobj_keep$subtypes014, levels = sub107order)
scPioneer::DimPlot_idx(Mobj_keep, group.by = "subtypes014", reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_subtype.pdf"),width = 18, height = 6)

pheatmap(colPercents(table(Mobj_keep$subtypes014, Mobj_keep$refine_ct03))[1:107,], cluster_rows = F, cluster_cols = F)
ggsave(filename = paste0(outdir,"/Mobj_subtypeXcelltype.pdf"),width = 8, height = 21)

tmpprop=colPercents(table(Mobj_keep$MajorClass02,Mobj_keep$Location_short))[1:6,]
tmpprop=tmpprop/rowMeans(tmpprop)
tmpprop=tmpprop %>% melt()
tmpprop=left_join(tmpprop,
                  data.frame(unique(Mobj_keep@meta.data[,c("Location_short","System")])),
                  by=c("Var2"="Location_short"))
tmpprop_filter=rbind(
  tmpprop %>% group_by(System,Var1) %>% top_n(3, value),
  tmpprop %>% group_by(System,Var1) %>% top_n(3, -value)
) %>% arrange(System, Var1, desc(value))
colnames(tmpprop_filter)=c("Cell_class","Location","FC","System")
tmpprop_filter$logFC=log2(tmpprop_filter$FC)

for(x in levels(Mobj_keep$MajorClass02)){
lollipop_plot(
  tmpprop_filter[tmpprop_filter$Cell_class %in% x,] %>% unique(),
  "Location","logFC",group = "System",
)#+theme_minimal_vgrid()+NoLegend()
ggsave(filename = paste0(outdir,"/Mobj_CellClass_LocEnrich_",gsub("\\.","",gsub("\\s+\\S+","",x)),".pdf"),width = 4.5, height = 6.5)
}

newColor=hue_pal()(length(Mobj_keep$System %>% levels()))
newColor=berryFunctions::addAlpha(newColor,0.5)
show_col(newColor)
names(newColor)=Mobj_keep$System %>% levels()

tmpprop00=colPercents(table(Mobj_keep$MajorClass02,Mobj_keep$System))[1:6,]
tmpprop00=tmpprop00/rowMeans(tmpprop00)
tmpprop00=tmpprop00 %>% melt()
tmpprop_filter00=tmpprop00 %>% arrange(Var1, desc(value))
colnames(tmpprop_filter00)=c("Cell_class","System","FC")
tmpprop_filter00$logFC=log2(tmpprop_filter00$FC+1)
for(x in levels(Mobj_keep$MajorClass02)){
  lollipop_plot(
    tmpprop_filter00[tmpprop_filter00$Cell_class %in% x,] %>% unique(),
    "System","logFC",
    dot_color = newColor
  )+theme(axis.line.y = element_line(colour = "grey"))+geom_vline(xintercept = 0, color = "lightgrey",linewidth=0.6)+geom_vline(xintercept = 1) #+NoAxes()
  #+theme_minimal_vgrid()+NoLegend()
  ggsave(filename = paste0(outdir,"/Mobj_SysEnrich_",gsub("\\.","",gsub("\\s+\\S+","",x)),".OnlySys.pdf"),width = 2.75, height = 2.75)
}

for(x in 1:length(target_var)){
  roveRun(Mobj_keep, Group=target_var[x], CT="MajorClass02", outdir=paste0(outdir,"/",target_var_save[x],"_CellClass"), flip = FALSE, mywidth=3+0.3*length(unique(Mobj_keep@meta.data[,target_var[x]])), nCT=10)
}
for(x in 1:length(target_var)){
  roveRun(Mobj_keep, Group=target_var[x], CT="refine_ct03", outdir=paste0(outdir,"/",target_var_save[x],"_Celltype"), flip = FALSE, mywidth=3+0.3*length(unique(Mobj_keep@meta.data[,target_var[x]])), nCT=35)
}
for(x in 1:length(target_var)){
  roveRun(Mobj_keep, Group=target_var[x], CT="subtypes014", outdir=paste0(outdir,"/",target_var_save[x],"_subtype"), flip = FALSE, mywidth=3+0.3*length(unique(Mobj_keep@meta.data[,target_var[x]])), nCT=107)
}
for(x in 1:length(target_var)){
  indf=Mobj_keep@meta.data[,c(target_var[x],"MajorClass02")]
  indf[[target_var_save[x]]]=indf[[target_var[x]]]
  indf=indf[,c(target_var_save[x],"MajorClass02")]
  scPioneer::plot_fraction(indf, color.use = hue_pal()(6))
  ggsave(filename = paste0(outdir,"/",target_var_save[x],"_CellClass_barplot.pdf"),width = 5.5, height = 2+0.2*length(unique(Mobj_keep@meta.data[,target_var[x]])))
}
for(x in 1:length(target_var)){
  indf=Mobj_keep@meta.data[,c(target_var[x],"refine_ct03")]
  indf[[target_var_save[x]]]=indf[[target_var[x]]]
  indf=indf[,c(target_var_save[x],"refine_ct03")]
  scPioneer::plot_fraction(indf, color.use = hue_pal()(32))
  ggsave(filename = paste0(outdir,"/",target_var_save[x],"_Celltype_barplot.pdf"),width = 7.5, height = 2+0.2*length(unique(Mobj_keep@meta.data[,target_var[x]])))
}

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
scPioneer::DotPlot2(Mobj_keep,features = Marker_Lev2, group.by = "refine_ct03")
ggsave(filename = paste0(outdir,"/Mobj_Marker_Lev2_DotPlot.pdf"),width = 11, height = 6)

rdf=rowPercents(table(Mobj_keep$Location_short, Mobj_keep$subtypes014))[1:51,1:107]
rdf=rdf[, colnames(rdf) %!in% levels(Mobj_keep$subtypes014)[32:61]]
rdf=rdf[!is.na(rowSums(rdf)),]
rdf_occ=(rdf>=1.5)+0
dfin=data.frame()
for(x in rownames(rdf_occ)){
  xxx=rdf_occ[x,]
  xxx=xxx[xxx>0]
  tmpdf=data.frame(cluster=rep(x, length(xxx)),gene=names(xxx))  ###gene here represents celltype
  dfin=rbind(dfin,tmpdf)
}
dfin_x=dfin
dfin_x$group=rep("group1",dim(dfin_x)[1])
dfin_y=dfin
dfin_y$group=rep("group2",dim(dfin_y)[1])
dfin2=rbind(dfin_x,dfin_y)
JACdfin=CorTable(dfin2, JC = TRUE)
pheatmap(JACdfin, cluster_rows = F, cluster_cols = F, display_numbers = JACdfin, name = "JacDis.")

### 20250117
outdir="./AnaRes20250117"
dir.create(outdir)

CT32_color=sample(scPioneer::scPalette2(length(levels(Mobj_keep$refine_ct03))))
names(CT32_color)=levels(Mobj_keep$refine_ct03)
show_col(CT32_color)
CT6_color=sample(scPioneer::scPalette2(length(levels(Mobj_keep$MajorClass02))))
names(CT6_color)=levels(Mobj_keep$MajorClass02)
CT6_color["Germ Cells"]="#DACA00"
show_col(CT6_color)
DimPlot(Mobj_keep, group.by = "MajorClass02", reduction = "umap", cols=CT6_color, pt.size = 1, ncol = 1, repel =T, label = FALSE)+ThemeBox+xlab("")+ylab("")+ggtitle("")
ggsave(filename = paste0(outdir,"/Mobj_cellclass.pdf"),width = 5.5, height = 4.25)
scPioneer::DimPlot_idx(Mobj_keep, group.by = "refine_ct03", reduction = "umap",  cols=as.vector(CT32_color), pt.size = 1, ncol = 1, repel =T)+ThemeBox+xlab("")+ylab("")
ggsave(filename = paste0(outdir,"/Mobj_celltype.pdf"),width = 8, height = 5)

Mobj_keep$subtypes015=Mobj_keep$subtypes014
Mobj_keep$subtypes015=gsub("lipoFB01_APOE","FB14_APOE",Mobj_keep$subtypes015)
Mobj_keep$subtypes015=gsub("myoFB01_ELN","FB15_ELN",Mobj_keep$subtypes015)
Mobj_keep$subtypes015=gsub("FB14_WT1","FBL01_WT1",Mobj_keep$subtypes015)
Mobj_keep$subtypes015=factor(Mobj_keep$subtypes015, levels =
                               c(levels(Mobj_keep$subtypes014)[1:87],
                                 "FB14_APOE","FB15_ELN","FBL01_WT1",
                                 levels(Mobj_keep$subtypes014)[91:107]))

CT107_color=sample(scPioneer::scPalette2(length(levels(Mobj_keep$subtypes015))))
names(CT107_color)=levels(Mobj_keep$subtypes015)
scPioneer::DimPlot_idx(Mobj_keep, group.by = "subtypes015", cols=as.vector(CT107_color), reduction = "umap", pt.size = 1, ncol = 1, label.idx = F)+
  ThemeBox+xlab("")+ylab("")#+scale_color_manual(values = CT107_color)
ggsave(filename = paste0(outdir,"/Mobj_subtype.pdf"),width = 18, height = 6)

#xxx=DimPlot_idx(dData_Imm, group.by = "CT89", prefix.index = 58:89, cols = as.vector(CT89_color[levels(dData_Imm$CT89)]), raster = T, raster.dpi=c(1024*0.7,1024*0.7), repel =T)+NoLegend()
iCT=levels(Mobj_keep$subtypes015)
nCT=length(iCT)
mycolortmp=CT107_color
show_col(mycolortmp)
names(mycolortmp)=iCT

ImmuData=subset(Mobj_keep, cellid %in% c(colnames(reSeuList$seu_immu),
                  colnames(subset(Mobj_keep, subtypes015 %in% c("Mast")))))
ImmuData=fastRefine(ImmuData, reduct = "pca_harmony", dims = 1:30, resol = 0.3, alg = 2, nnei = 150, mindist = 0.5, spread = 1.25)
ImmuData=subset(ImmuData, seurat_clusters %!in% c(7))
xxxumap=ImmuData@reductions$umap@cell.embeddings
yyyumap=reSeuList$seu_immu@reductions$umap@cell.embeddings
for(x in rownames(xxxumap)){
  if(x %in% rownames(yyyumap)){
  xxxumap[x,"umap_1"]=yyyumap[x,"umap_1"]
  xxxumap[x,"umap_2"]=yyyumap[x,"umap_2"]
  }else{
    xxxumap[x,"umap_1"]=xxxumap[x,"umap_1"]+3.5
    xxxumap[x,"umap_2"]=xxxumap[x,"umap_2"]-1.25
  }
}
ImmuData[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(xxxumap), key = "umap_", assay = DefaultAssay(ImmuData))
ImmuData$subtypes015=droplevels(ImmuData$subtypes015)
scPioneer::DimPlot_idx(ImmuData, group.by = "subtypes015", prefix.index = 1:31,
            cols = as.vector(mycolortmp[levels(ImmuData$subtypes015)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
ggsave(filename = paste0(outdir,"/Immu_subtype.pdf"),width = 4.5, height = 4.5)

seu_epi_curat$forCombine02=seu_epi_curat$forCombine01
seu_epi_curat=modifyAno(seu_epi_curat, anoCol="forCombine02", oldID = "^CNTN4_Cells", newID = "CNTN4_NEL_Cells",
                    newOrder = gsub("^CNTN4_Cells","CNTN4_NEL_Cells",levels(seu_epi_curat$forCombine02)))
seu_epi_curat=subset(seu_epi_curat, forCombine02 %!in% c("Glia_Cells"))
seu_epi_curat$forCombine02=droplevels(seu_epi_curat$forCombine02)
scPioneer::DimPlot_idx(seu_epi_curat, group.by = "forCombine02", prefix.index = 32:61,
            cols = as.vector(mycolortmp[levels(seu_epi_curat$forCombine02)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
ggsave(filename = paste0(outdir,"/Epi_subtype.pdf"),width = 4.5, height = 4.5)

scPioneer::DimPlot_idx(reSeuList$seu_endo, group.by = "subtypes00", prefix.index = 62:72,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_endo$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.2,1024*0.2), repel =T)+NoLegend()+NoAxes()
ggsave(filename = paste0(outdir,"/Endo_subtype.pdf"),width = 4.5, height = 4.5)

reSeuList$seu_strom$subtypes01=reSeuList$seu_strom$subtypes00
reSeuList$seu_strom$subtypes01=gsub("lipoFB01_APOE","FB14_APOE",reSeuList$seu_strom$subtypes01)
reSeuList$seu_strom$subtypes01=gsub("MyoFB01_ELN","FB15_ELN",reSeuList$seu_strom$subtypes01)
reSeuList$seu_strom$subtypes01=gsub("MESO01_WT1","FBL01_WT1",reSeuList$seu_strom$subtypes01)
reSeuList$seu_strom$subtypes01=factor(reSeuList$seu_strom$subtypes01, levels =
                               c(levels(reSeuList$seu_strom$subtypes00)[1:15],
                                 "FB14_APOE","FB15_ELN","FBL01_WT1",
                                 levels(reSeuList$seu_strom$subtypes00)[19:25]))
scPioneer::DimPlot_idx(reSeuList$seu_strom, group.by = "subtypes01", prefix.index = 73:97,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_strom$subtypes01)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
allfb=subset(reSeuList$seu_strom, subtypes01 %!in% c("SMC01_MCAM","SMC02_RERGL","SMC03_KALRN",
                                                    "SMC04_TAGLN","SMC05_RGS6","SMC06_AOC3","PER01_PDGFRB"))
allfb$subtypes01=droplevels(allfb$subtypes01)
table(allfb$subtypes01)

allfb$Organ_bk=allfb$Organ
allfb$Organ=allfb@meta.data[,"Location_short"]
for(i in 1:length(refine_list)){
  for(j in 1:length(refine_list[[i]])){
    allfb=modifyAno(allfb, anoCol = "Organ", oldID = paste0("^",refine_list[[i]][j],"$"), newID = names(refine_list[i]))
  }
}
allfb$Organ=factor(allfb$Organ, levels = rev(names(refine_list)))
allfb$Organ=droplevels(allfb$Organ)

allfb@reductions$umap_bak=allfb@reductions$umap
obsm_ct=allfb@meta.data
obsm=subdataList$strom_fb@reductions$umap@cell.embeddings %>% as.data.frame()
obsm=obsm[rownames(allfb@reductions$umap_bak@cell.embeddings),] #%>% as.matrix()
obsm$cellid=rownames(obsm)
obsm_ct=left_join(obsm_ct, obsm, by="cellid")
#allfb[['umap']] <- CreateDimReducObject(embeddings = obsm %>% as.matrix(), assay = "RNA", key = "umap_")
ggplot(obsm_ct)+geom_point(aes(umap_1, umap_2, colour = subtypes01), size=0.1)+
  scale_color_manual(values = scPioneer::scPalette2(length(unique(obsm_ct$subtypes01))))+
  xlim(c(-20,20)) + ylim(c(-15,10)) + theme_minimal() + theme(panel.grid = element_blank()) +NoAxes()
ggplot()+stat_density_2d(data=obsm_ct,
                         aes(umap_1, umap_2, color=as.factor("..level..")), breaks=c(0.0002,100), h = c(3, 3), n=100)+
  scale_color_manual(values = c("black",rep("NA",100),"NA"))+
  stat_density_2d(data=obsm_ct,
                  aes_string(x = "umap_1", y = "umap_2",  alpha = "..level..", fill = "subtypes01"),
                  linewidth = 0, geom = "polygon",  breaks=c(0.01,50), n = 200, h = c(0.9, 0.9)) + theme_cowplot() +
  scale_fill_manual(values = scPioneer::scPalette2(length(unique(obsm_ct$subtypes01))))+
  NoAxes() + xlim(c(-20,20)) + ylim(c(-15,10))

roveRun(allfb, Group="Organ", CT="subtypes01", outdir=paste0("./data/LocationFBtypes"), flip = FALSE, mywidth=15, nCT=20)
DimPlot(allfb, reduction = "umap")

scPioneer::DimPlot_idx(allfb, group.by = "subtypes01", #prefix.index = 73:97,
                       raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()

ggsave(filename = paste0(outdir,"/Stroma_subtype.pdf"),width = 4.5, height = 4.5)

scPioneer::DimPlot_idx(reSeuList$seu_germ, group.by = "subtypes00", prefix.index = 98:106,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_germ$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.2,1024*0.2), repel =T)+NoLegend()+NoAxes()
ggsave(filename = paste0(outdir,"/Germ_subtype.pdf"),width = 4.5, height = 4.5)

Mobj_tmp=subset(Mobj_keep, subtypes015 %!in% c("Prolif."))
Mobj_tmp$subtypes015=droplevels(Mobj_tmp$subtypes015)
xxx=table(Mobj_tmp$System, Mobj_tmp$subtypes015)
yyy=colPercents(xxx)[1:dim(xxx)[1],]
yyy=melt(yyy)
colnames(yyy)=c("group","CT","proportion")
yyy$CT=factor(yyy$CT, levels = rev(levels(Mobj_tmp$subtypes015)))
yyy$group=factor(yyy$group, levels = levels(Mobj_tmp$System))

zzz1=yyy
zzz1$group=if_else(zzz1$group %in% c("DS"),"DS", "other")
addDF1=aggregate(proportion ~ group + CT, data = zzz1, FUN = sum, na.rm = TRUE)
addDF1=addDF1[addDF1$group %in% c("DS"),]
zzz2=yyy
zzz2$group=if_else(zzz2$group %in% c("DS","RS"),"RS", "other")
addDF2=aggregate(proportion ~ group + CT, data = zzz2, FUN = sum, na.rm = TRUE)
addDF2=addDF2[addDF2$group %in% c("RS"),]
zzz3=yyy
zzz3$group=if_else(zzz3$group %in% c("DS","RS","US"),"US", "other")
addDF3=aggregate(proportion ~ group + CT, data = zzz3, FUN = sum, na.rm = TRUE)
addDF3=addDF3[addDF3$group %in% c("US"),]
zzz4=yyy
zzz4$group=if_else(zzz4$group %in% c("DS","RS","US","RPS"),"RPS", "other")
addDF4=aggregate(proportion ~ group + CT, data = zzz4, FUN = sum, na.rm = TRUE)
addDF4=addDF4[addDF4$group %in% c("RPS"),]
zzz5=yyy
zzz5$group=if_else(zzz5$group %in% c("DS","RS","US","RPS","CS"),"CS", "other")
addDF5=aggregate(proportion ~ group + CT, data = zzz5, FUN = sum, na.rm = TRUE)
addDF5=addDF5[addDF5$group %in% c("CS"),]

yyy_plot=ggplot(yyy)+geom_bar(aes(CT, proportion, group=group, fill=group), width=0.5, stat = "identity", position = position_stack())+
  geom_point(data=addDF1, aes(CT, 100-proportion), shape=124, size=1, color="black")+   #shape=95: -; shape=124: |
  geom_line(data=addDF1, aes(CT, 100-proportion, group=1),  stat = "identity", position = "identity", linetype="solid", color="black")+
  geom_point(data=addDF2, aes(CT, 100-proportion), shape=124, size=1, color="darkgray")+   #shape=95: -; shape=124: |
  geom_line(data=addDF2, aes(CT, 100-proportion, group=1),  stat = "identity", position = "identity", linetype="solid", color="darkgray")+
  geom_point(data=addDF3, aes(CT, 100-proportion), shape=124, size=1, color="gray")+   #shape=95: -; shape=124: |
  geom_line(data=addDF3, aes(CT, 100-proportion, group=1),  stat = "identity", position = "identity", linetype="solid", color="gray")+
  geom_point(data=addDF4, aes(CT, 100-proportion), shape=124, size=1, color="lightgray")+   #shape=95: -; shape=124: |
  geom_line(data=addDF4, aes(CT, 100-proportion, group=1),  stat = "identity", position = "identity", linetype="solid", color="lightgray")+
  geom_point(data=addDF5, aes(CT, 100-proportion), shape=124, size=1, color="white")+   #shape=95: -; shape=124: |
  geom_line(data=addDF5, aes(CT, 100-proportion, group=1),  stat = "identity", position = "identity", linetype="solid", color="white")+
  scale_fill_manual(values = hue_pal()(6))+theme_minimal()+theme(axis.text = element_blank())+ #theme(axis.text = element_text(size=14, colour = "black"))+
  xlab("")+ylab("")+coord_flip()+NoLegend()
runFig(paste0(outdir,"/barplot.6Sys105CT_NoProlif_NoLegend"), yyy_plot, 3, 18)

ncols=1
xxxplot=scPioneer::DimPlot_idx(Mobj_tmp, group.by = "subtypes015", idx.sep = "  ", prefix.index = NULL,
                    cols=as.vector(mycolortmp[names(mycolortmp) %!in% c("Prolif.")]))+
  guides(col=guide_legend(ncol = ncols, byrow = FALSE, title=NULL))
xxx=get_legend(xxxplot)
npoints=sum(grepl("points",xxx[["grobs"]][[1]]$grobs)==TRUE)
lll=xxx[["grobs"]][[1]][["grobs"]]
widthV=xxx[["grobs"]][[1]]$widths
widthV_max=max(widthV)
for(ncolidx in 1:ncols){
  widthV[ncolidx*4]=widthV_max
}
xxx[["grobs"]][[1]]$widths=widthV
xxx$widths[3]=sum(widthV)
point_idx=0
yyy=xxx
for(li in 1:length(lll)){
  if(grepl("points",lll[[li]])){
    lll[[li]][["size"]]=unit(0.65, "cm")
    lll[[li]][["x"]]=unit(1, "native") #3.5*lll[[li]][["x"]]#unit(0.65, "cm")
  }
  if(grepl("rect",lll[[li]])){
    lll[[li]]$gp$fill="white"
  }
  if(grepl("gTree",lll[[li]])){
    point_idx=point_idx+1
    lll[[li]]$vp$x=unit(0, "npc")

    tmplab=lll[[li]]$children[[1]]$children[[1]]$label
    tmplab1=strsplit(tmplab, split = "  ")[[1]][1]
    tmplab2=strsplit(tmplab, split = "  ")[[1]][2]
    if(as.numeric(tmplab1)<10){
      lll[[li]]$children[[1]]$children[[1]]$label=paste0(rep("", 1), rep("\t-", 1),tmplab2)
    }else{
      lll[[li]]$children[[1]]$children[[1]]$label=paste0(rep("\t-",1),tmplab2)
    }
  }
}
xxx[["grobs"]][[1]][["grobs"]]=lll
Aplot=plot_grid(NULL,
                plot_grid(NULL,xxx,nrow = 3, rel_heights = c(0.01,0.98,0.01)),
                NULL,
                nrow = 1, rel_widths = c(0.12,0.87,0.01))
## second layer
npoints=sum(grepl("points",yyy[["grobs"]][[1]]$grobs)==TRUE)
lll=yyy[["grobs"]][[1]][["grobs"]]
point_idx=0
for(li in 1:length(lll)){
  if(grepl("points",lll[[li]])){
    lll[[li]][["size"]]=unit(0, "cm")
    lll[[li]][["x"]]=unit(1, "native")
  }
  if(grepl("rect",lll[[li]])){
    lll[[li]]$gp$fill="white"
  }
  if(grepl("gTree",lll[[li]])){
    point_idx=point_idx+1

    lll[[li]]$vp$x=unit(-0.065, "npc")
    tmplab=lll[[li]]$children[[1]]$children[[1]]$label
    tmplab1=strsplit(tmplab, split = "  ")[[1]][1]
    tmplab2=strsplit(tmplab, split = "  ")[[1]][2]
    lll[[li]]$children[[1]]$children[[1]]$label=
      if_else(as.numeric(tmplab1)<10, paste0(rep("\t", 1),tmplab1), tmplab1) #tmplab1 #paste0(rep("", 1),tmplab1)

    tmpcolor=lll[[li-2*npoints+point_idx]]$gp$col
    hcl=farver::decode_colour(tmpcolor, "rgb", "hcl")
    label_col <- ifelse(hcl[, "l"] > 50, "black", "lightgrey")
    lll[[li]]$children[[1]]$children[[1]]$gp$col=as.vector(label_col)
  }
}
yyy[["grobs"]][[1]][["grobs"]]=lll
Bplot=plot_grid(NULL,
                plot_grid(NULL,yyy,nrow = 3, rel_heights = c(0.01,0.98,0.01)),
                NULL,
                nrow = 1, rel_widths = c(0.12,0.87,0.01))
ppp=Aplot+annotation_custom(grob = ggplotGrob(Bplot))
runFig(paste0(outdir, "/All106.6Sys106subCT_NoProlif.CTlegend"), ppp, 4.25, 25)

Location_color=hue_pal()(length(levels(Mobj_keep$Location_short)))
names(Location_color)=levels(Mobj_keep$Location_short)

yyy=colPercents(table(Mobj_keep$MajorClass02,Mobj_keep$Location_short))[1:4,] %>% as.data.frame()
yyy_order=rownames(yyy)
yyy$Celltype=rownames(yyy)
yyy_long=melt(yyy,id.vars = "Celltype")
colnames(yyy_long)=c("Celltype","Location","Rate")
yyy_long=yyy_long %>% arrange(Celltype, desc(Rate))
yyy_long=yyy_long %>% group_by(Celltype) %>% top_n(n=15, Rate)

data=yyy_long %>% as.data.frame()

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*length(yyy_order), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Celltype <- rep(yyy_order, each=empty_bar)
data <- rbind(data, to_add)
data$Celltype=factor(data$Celltype,yyy_order)
data=data[order(data$Celltype),]
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

base_data <- data %>%
  group_by(Celltype) %>%
  summarize(start=min(id), end=max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))
# prepare a data frame for grid (scales)

# Make the plot
ggplot(data, aes(x=as.factor(id), y=Rate, fill=Location)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar

  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.

  # Add text showing the value of each 100/75/50/25 lines

  geom_bar(stat="identity", alpha=0.5) +
  scale_fill_manual(values = Location_color)+
  ylim(-120,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=Rate+10, label=Location, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -25, label=Celltype), hjust=c(0.5,0.5,0.5,0.5),  angle=c(-40,50,-40,50), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(filename = paste0(outdir,"/location_subtype_top15.pdf"),width = 5.5, height = 5.5, scale = 1.1)

idx=1:length(CT107_color)
names(idx)=names(CT107_color)

CT6toCT106=unique(Mobj_keep@meta.data[,c("MajorClass02","subtypes015")])
targetCT="Stromal Cells"
top6=6

refdata=yyy_long[yyy_long$Celltype==targetCT,]
refdata_sub=CT6toCT106[CT6toCT106$MajorClass02==targetCT,]
selectedData=colPercents(table(Mobj_keep$subtypes015,Mobj_keep$Location_short))[levels(Mobj_keep$subtypes015)[levels(Mobj_keep$subtypes015) %in% refdata_sub$subtypes015],
                                                                                refdata$Location]
selectedCT=selectedData %>% rowSums() %>% sort() %>% tail(n=top6) %>% names()
selectedData=selectedData[rev(selectedCT),] %>% as.data.frame()
selectedData$idx=as.character(idx[rownames(selectedData)])
zzz= selectedData %>% melt(id.vars = "idx")

#tmpPlot+coord_polar(theta = "y", direction = -1)+theme_nothing()

zzz_add=data.frame(idx=lapply(unique(zzz$idx),FUN = function(x){rep(x,nrow(data))}) %>% unlist(),
                   variable=rep(data$Location,top6),
                   value=rep(data$Rate,top6),
                   ct=rep(data$Celltype,top6),
                   id=rep(data$id, top6))
zzz_add$value=0
for(x in 1:nrow(zzz_add)){
  zzztmp=zzz_add[x,]
  if(zzztmp$ct==targetCT &
     zzztmp$idx %in% zzz$idx &
     zzztmp$variable %in% zzz$variable){
    zzz_add[x, "value"]=zzz[zzz$idx == zzztmp$idx & zzz$variable == as.character(zzztmp$variable),"value"]
  }
}
zzz_add$idx=factor(zzz_add$idx, levels = rev(as.character(sort(as.numeric(unique(zzz_add$idx))))))
zzz_add$id=factor(zzz_add$id, levels = as.character(sort(as.numeric(unique(zzz_add$id)))))

ppp=ggplot(data = zzz_add, aes(x=id, y=idx, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank()) + ###, axis.text.y=element_blank(),
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+NoLegend()
ppp+coord_polar()
ggsave(filename = paste0(outdir,"/location_subtype_top15loc_top6subtype_",targetCT,".pdf"),width = 5.5, height = 5.5, scale = 1.1)
