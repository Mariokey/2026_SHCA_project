# ============================================================================
# Organ-wise ligand-receptor inference, network centrality and visualisation
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Cell-cell interactions"
# Source     : extracted from Analysis240525.r (L7863-8150, L11050-11116)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### CCI subtypes 20241011
#subsysList$DS/subtypes012
#subdataList$strom_fb/FB_SubClass02
#subdataList$myl_mp/MYL_SubClass04_bak

CCIdir_img="./data/R_data/CCI_subDS/"
dir.create(CCIdir_img, recursive = T)
Mobj_DS=subset(Mobj, Organ %in% c("Liver","BS","UDT","Intestine")) ###64663 cells
intersect_cb=intersect(Mobj_DS$cellid_ori %>% as.vector(),
                       c(subsysList$DS$cellid %>% as.vector(),
                         subdataList$strom_fb$cellid %>% as.vector(),
                         subdataList$myl_mp$cellid %>% as.vector())
                       )
Mobj_DS_oi=subset(Mobj_DS, cellid_ori %in% intersect_cb)

ds_episub=subset(subsysList$DS, cellid %in% as.vector(Mobj_DS_oi$cellid_ori))
ds_mpsub=subset(subdataList$myl_mp, cellid %in% as.vector(Mobj_DS_oi$cellid_ori))
ds_fbsub=subset(subdataList$strom_fb, cellid %in% as.vector(Mobj_DS_oi$cellid_ori))

Nsubobjlist=list()
Nsubobjlist["ds_episub"]=list(c("subtypes012"))
Nsubobjlist["ds_mpsub"]=list(c("MYL_SubClass04_bak"))
Nsubobjlist["ds_fbsub"]=list(c("FB_SubClass02"))

Mobj_DS_oi$subtypes012=Mobj_DS_oi$subtypes011
Mobj_DS_oi$subtypes012=as.character(Mobj_DS_oi$subtypes012)
Mobj_DS_oi=addNewMetaColumn(Mobj_DS_oi, Nsubobjlist, nameV=NULL, forceAdd=TRUE, forceAddname="subtypes012")

Mobj_DS_oi=subset(Mobj_DS_oi, subtypes012 %!in% c("rmv_BS","rmv_Intestine","rmv_UDT"))

suborder012=c(levels(ds_episub$subtypes012),
              levels(ds_mpsub$MYL_SubClass04_bak),
              levels(ds_fbsub$FB_SubClass02))
suborder012=suborder012[!grepl("rmv_",suborder012)]
Mobj_DS_oi$subtypes012=factor(Mobj_DS_oi$subtypes012, levels=suborder012)
Mobj_DS_oi$subtypes012=droplevels(Mobj_DS_oi$subtypes012)
Mobj_DS_oi$Organ=droplevels(Mobj_DS_oi$Organ)

xxxdf=colPercents(table(Mobj_DS_oi$subtypes012, Mobj_DS_oi$Organ))

for(xi in 1:length(CCI_neg_list)){
  x=names(CCI_neg_list[xi])
  rmv_ct=CCI_neg_list[[xi]]

  Mobj_tmp=subset(Mobj_DS_oi, Organ %in% c(x))
  Mobj_tmp=subset(Mobj_tmp, subtypes012 %!in% c(rmv_ct))
  Mobj_tmp$Organ=droplevels(Mobj_tmp$Organ)
  Mobj_tmp$subtypes012=droplevels(Mobj_tmp$subtypes012)
  print(x)
  print(table(Mobj_tmp$subtypes012, Mobj_tmp$Organ))

  Mobj_tmp$CCItype=Mobj_tmp$subtypes012
  runCCI_perData(Mobj_tmp,
                 spatialRun = FALSE, #inlist = st_preCfg,
                 group = c(x), groupby = "Organ", indir = CCIdir_img, drawOnly=FALSE)
}

CCIrds=list.files(CCIdir_img)
CCIrds=CCIrds[grepl(".rds",CCIrds)]
allCCIlr_DF=data.frame()
for(xi in CCIrds){
  xi_data=readRDS(paste0(CCIdir_img,"/",xi))
  xi=gsub(".CCI.rds","",xi)

  idx=levels(xi_data@idents)
  bulb=netVisual_bubble(xi_data, sources.use = c(1:length(idx)), targets.use = c(1:length(idx)),
                        angle.x = 45, return.data=T)
  tmpDF=bulb$communication
  tmpDF$tissue=xi

  allCCIlr_DF=rbind(allCCIlr_DF, tmpDF)
}
write.xlsx(allCCIlr_DF, file = paste0(CCIdir_img,"/","DS_subOI_CCIlr_DF.across_tissue.xlsx"))

################################################################################

###
outdir="./data/R_data/AnaRes20240907/"
dir.create(outdir)
roveRun(inobj, Group="Location_short", CT="MYL_SubClass04",
        outdir=outdir, flip = FALSE, mywidth=8, nCT=15)

xxxn$Group=droplevels(xxxn$Group)
roveRun(xxxn, Group="Group1", CT="ident",         outdir="./rove.pdf", flip = FALSE, mywidth=8, nCT=5)

roveRun(subdataList$strom_fb, Group="Location_short", CT="FB_SubClass02",
        outdir=paste0(outdir,"/FB_"), flip = FALSE, mywidth=8, nCT=15)

###
tmpdir="DATA_DIR"
infiles=list.files(tmpdir)
infiles=infiles[grepl("data_",infiles)]
refAnnoDF=read.xlsx(paste0(tmpdir, "ori_arrange.xlsx" ))
refAnnoDF$join_col=paste0(refAnnoDF$Tumor_type_original,"_",refAnnoDF$PosID)

xxxMark="EC04"
xxxMark="MP08"
xxxMark="per01"
xxxMark="per02"

inTest=infiles[grepl(xxxMark,infiles)]

AllExpDF=data.frame()
for(xlsx in inTest){
ncol_names=c("Name", "total", "g1", "g2", "g12")
sheet_names <- openxlsx::getSheetNames(paste0(tmpdir, '/', xlsx))
for(i in sheet_names){
iii=gsub("^.+\\s+","",i)
tmpdf=openxlsx::read.xlsx(paste0(tmpdir, '/', xlsx), sheet = i)
colnames(tmpdf)=ncol_names

tmpdf$label=lapply(1:nrow(tmpdf),function(i){unlist(strsplit(tmpdf$Name[i],"\\(|\\)|, "))[3]}) %>% unlist()
splitVs=lapply(1:nrow(tmpdf),function(i){unlist(strsplit(tmpdf$Name[i],"\\(|\\)|, "))[2]}) %>% unlist()
for(x in 1:length(splitVs)){splitVs[x]=if_else(nchar(splitVs[x])==1, paste0("0",splitVs[x]), splitVs[x])}
tmpdf$PosID=paste0(tmpdf$label, splitVs)
tmpdf$g1_pct=100*(tmpdf$g1+tmpdf$g12)/tmpdf$total
tmpdf$g2_pct=100*(tmpdf$g2+tmpdf$g12)/tmpdf$total
tmpdf$g12_pct=100*tmpdf$g12/tmpdf$total
tmpdf$g12_g1_pct=100*tmpdf$g12/(tmpdf$g1+tmpdf$g12)
tmpdf$g12_g2_pct=100*tmpdf$g12/(tmpdf$g2+tmpdf$g12)
tmpdf$Tumor_type_original00=iii
#  }
#}
AllExpDF=rbind(AllExpDF, tmpdf)
}
}
AllExpDF$join_col=paste0(AllExpDF$Tumor_type_original00,"_",AllExpDF$PosID)
AllExpDF_join=left_join(AllExpDF, refAnnoDF, by="join_col")

write.xlsx(AllExpDF_join, paste0(tmpdir, "/",xxxMark,"_allcombined.xlsx"))

myvectors=list()
myvectors[["TMA1"]]=c("TMA1","TMA")
myvectors[["TMA2"]]=c("TMA2","TMA")
myvectors[["TMA12"]]=c("TMA1","TMA2","TMA")

for(myvectors_idx in 1:length(myvectors)){
myvectors_names=names(myvectors[myvectors_idx])
plotdata=AllExpDF_join[AllExpDF_join$Region %in% c("nLUN","tLUN","bmBRN") &
                       #AllExpDF_join$Tumor_type_original %in% c("TMA1","TMA2","TMA") &
                       #AllExpDF_join$Tumor_type_original %in% c("TMA1","TMA") &
                       AllExpDF_join$Tumor_type_original %in% myvectors[[myvectors_names]] &
            AllExpDF_join$Tumor_type %in% c("ADC",
                                            "ADC_nLUN"),]
plotdata$Region=factor(plotdata$Region, levels = c("nLUN","tLUN","bmBRN"))

draw_plotlist=list()
for(item in c("total","g1","g2","g12","g1_pct","g2_pct","g12_pct","g12_g1_pct","g12_g2_pct")){

plotdata$exp=plotdata[,item]
tmax=max(plotdata$exp)
tmin=min(plotdata$exp)
plotdata$group=plotdata$Region

plot1=plotpubr(plotdata, gfacet = NULL, cmp=list(c("nLUN", "tLUN"),c("nLUN","bmBRN"),c("tLUN","bmBRN")))+
NoLegend()+xlab("")+ylab(item)+myfont+ylim(tmin,tmax*(1.5))
draw_plotlist[[item]]=plot1

}
xxx=plot_grid(plotlist = draw_plotlist, ncol = 1, scale = 1, align=c("hv"), axis=c("tblr"))
runFig(paste0(tmpdir,"/", xxxMark,"_data_",myvectors_names), xxx, 3.25, 24)
}

plotpubr=function(indf, gfacet=NULL, cmp=NULL, ...){

  splitdf1=ggpubr::compare_means(data = indf, formula = exp~group, group.by = gfacet, method = "t.test")
  stat1=ggpubr::stat_compare_means(data = indf, aes(label = ..p.signif..), vjust=0, bracket.size=0.5, hide.ns = FALSE, comparisons = cmp, method = "t.test")

  hhsplitbox1=ggpubr::ggviolin(indf, x = "group", y = "exp", width = 0.8, scale="width",
                               color="group", fill="group", trim = TRUE,
                               add = c("boxplot"),
                               add.params=list(shape=21, color=c("black"), fill="white", size=0.3, width=0.3),
                               short.panel.labs = TRUE,  #facet.by = gfacet, nrow=1,
                               ggtheme =theme_minimal())

  hhsplitbox1=hhsplitbox1+stat1

  return(hhsplitbox1)

}

tmpcsv=read.csv("DATA_DIR")
head(tmpcsv)
ggplot(tmpcsv[
  tmpcsv$Samples<100 &
    tmpcsv$Samples>1 &
    is.na(tmpcsv$Note),])+
  geom_point(aes(Samples,ExeHours, color=PRJ_type))

xxx=subset(seu_epi_all, Location_short  %in% c("ITJM")) # %in% c("ITJM","ITIM"))
xxx$Donor=factor(xxx$Donor, levels = c("39","8","2"))
xxx_color=rev(scales::hue_pal()(3))
names(xxx_color)=c("2","8","39")
addgenes=c("SLC27A2","SLC27A4","KHK","REG3A","REG3G", "SLC2A2","SLC2A6", "NEUROD1", "CHGB","LCT","ALPI","SI","TMPRSS15")
addgenes=c("CD36","FABP1","FABP2","FABP3","FABP4","SLC2A2","SLC2A6", "NEUROD1", "CHGB","LCT","ALPI","SI","TMPRSS15") #"PIGR","PHGR1")
VlnPlot(subset(xxx, subtypes012 %in% c("Metallothionein_Cells","Enterocytes",
                                       "Goblet_Cells","BEST4_Epithelial_Cells",
                                       "Endocrine_Cells", "Tuft_Cells", "Epi_undef2")),
               features = addgenes, group.by = "subtypes012", split.by = "Donor", ncol = 4, pt.size = 0,stack = T)+ #+geom_boxplot(width=0.25)+NoLegend()+xlab("")
               scale_fill_manual(values=xxx_color)#+

tdata=subdataList$lym_tnk
addgenes=c("IL17A","IL22","CD3D","CD3E","CD4","CD40LG","CD8A","CD8B")
markergenes=gsub(".+_","",levels(tdata$LYM_SubClass00))
markergenes=markergenes[markergenes %!in% "hsp"]
llist=list()
for(x in c("ILC3_IL4I1")){
ppp=VlnPlot(subset(tdata, LYM_SubClass00 %in% c(x)),
        features = c(addgenes, markergenes[markergenes %!in% c("FOXP3")]), group.by = "Donor",  #group.by = "LYM_SubClass00",
        split.by = "Donor", pt.size = 0,stack = T)+NoLegend() #+geom_boxplot(width=0.25)+xlab("")
runFig(paste0("./MP_mn2/", x, ".Tfeatures_splitbyDonor_vln1024"), ppp, w = 12, h = 8)
llist[[x]]=ppp
}

addgenes=c("CD4","FOLR2","LYVE1","OLFM4")
VlnPlot(seuloc_ITs_Epi_trajec,
        features = addgenes, group.by = "seurat_clusters", split.by = "Donor", pt.size = 0,stack = T) #+geom_boxplot(width=0.25)+NoLegend()+xlab("")

CCIdir_img="./CCI_sub0310/"
dir.create(CCIdir_img)

  rmv_ct="Prolif."
  Mobj_tmp=subset(Mobj_keep, Organ %!in% c("Pancreas","Testis","VD","Ovary","FU"))
  Mobj_tmp=subset(Mobj_tmp, subtypes015 %!in% c(rmv_ct))

  Mobj_tmp$Organ=droplevels(Mobj_tmp$Organ)
  Mobj_tmp$subtypes015=droplevels(Mobj_tmp$subtypes015)

  Mobj_tmp$CCItype=Mobj_tmp$subtypes015
  runCCI_perData(Mobj_tmp,
                 spatialRun = FALSE, #inlist = st_preCfg,
                 group = levels(Mobj_tmp$Organ), groupby = "Organ", indir = CCIdir_img, drawOnly=FALSE)

CCIrds=list.files(CCIdir_img)
CCIrds=CCIrds[grepl(".rds",CCIrds)]
allCCIlr_DF=data.frame()
for(xi in CCIrds){
  xi_data=readRDS(paste0(CCIdir_img,"/",xi))
  xi=gsub(".CCI.rds","",xi)

  idx=levels(xi_data@idents)
  bulb=netVisual_bubble(xi_data, sources.use = c(1:length(idx)), targets.use = c(1:length(idx)),
                        angle.x = 45, return.data=T)
  tmpDF=bulb$communication
  tmpDF$tissue=xi

  allCCIlr_DF=rbind(allCCIlr_DF, tmpDF)
}
write.xlsx(allCCIlr_DF, file = paste0(CCIdir_img,"/","All_subOI_CCIlr_DF.across_tissue.xlsx"))

tPR03=xi_data
tPR03=netAnalysis_computeCentrality(tPR03)
CellChat::netAnalysis_signalingRole_scatter(tPR03)

CToi="Mph01_FABP4"  #"Mph02_RGS1"
rankNet(tPR03,mode = "single",sources.use = CToi, stacked = T, do.stat = TRUE)+
rankNet(tPR03,mode = "single",targets.use = CToi, stacked = T, do.stat = TRUE)

CToi="CNTN4_NEL_Cells"
netVisual_bubble(tPR03, sources.use = CToi)+coord_flip()
ggsave(paste0(CCIdir_img, "/", CToi, "_asS_bubble.pdf"), width = 24, height = 12)
netVisual_bubble(tPR03, targets.use = CToi)+coord_flip()
ggsave(paste0(CCIdir_img, "/", CToi, "_asT_bubble.pdf"), width = 35, height = 12)

source_path=c("THBS","GALECTIN","CXCL","LXA4","MIF")
target_path=c("SELPLG","THBS","ICAM","VWF","GALECTIN")

pdf(width = 4,height = 12,file = paste0(CCIdir_img,"/pathway_for_platelet.","as_target",".circle.05",'.pdf'))
par(mfrow = c(5,1), mar=c(0.3,0.3,1.5,0.3), xpd=TRUE)
for(ppp in target_path){
  netVisual_aggregate(tPR03, signaling = ppp,  targets.use = "Platelet", reduce = 1,
                      layout = "circle",  edge.width.max = 10, signaling.name = paste0(ppp))
}
dev.off()
