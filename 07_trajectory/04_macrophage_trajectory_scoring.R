# ============================================================================
# Macrophage trajectory scoring and integration of trajectory results into the atlas object
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Trajectory analysis"
# Source     : extracted from Analysis240525.r (L5846-5882, L8278-8338)
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

### NOTE DO NOT load scPioneer!!! ->boxplot will be replaced by this package!!!
plotCytoTRACE(xdata_results, phenotype = ptV,  outputDir = dir_img,
              emb = xdata@reductions$umap@cell.embeddings %>% as.data.frame()) #, gene = "Kit")

mn2.xdata_trajec_UDT.Int=runMonocle2(subset(xdata, Organ %in% c("UDT","Intestine")),
                             "MPsubtype00", mygroup = "Location_short", mode = c("disp"),
                            rmvIdent=NULL, savefix=paste0(dir_img,"/cds_mp_trajec_1133.","MPsubtype00_","UDT.Int","."))

mn2.xdata_trajec_UDT.Int=readRDS("./workspace/PRJ21_mOrganYang/R_data/MP_reScore241104/cds_mp_trajec_1133.MPsubtype00_UDT.Int.monocle_cds_plotlist.dispdeg.rds")

library(cowplot)

newMPdata=mn2.xdata_trajec_UDT.Int
newMPdata$disp$pseudotime$data$Pseudotime=max(newMPdata$disp$pseudotime$data$Pseudotime)-newMPdata$disp$pseudotime$data$Pseudotime
xdf=newMPdata$disp$pseudotime$data
xdf=xdf[!(xdf$MPsubtype00 %in% c("c3_FCGR3A+.inMo")),]
xdf$group="MoDM"
xdf[xdf$MPsubtype00 %in% c("c4_RGS1+.Mph","c5_FOLR2+.Mph"),]$group="FeDM"
ccc=names(table(newMPdata$disp$pseudotime$data$MPsubtype00))
ccc_color=scPioneer::scPalette2(length(ccc))
names(ccc_color)=ccc
ggplot()+
  geom_violin(data=xdf, aes(MPsubtype00, Pseudotime, fill=MPsubtype00), scale = "width", trim = FALSE, adjust=2)+
  scale_fill_manual(values = ccc_color)+
  geom_boxplot(data=xdf, aes(MPsubtype00, Pseudotime), width=0.15, fill="white", outlier.alpha = 0)+
  xlab("")+theme_cowplot()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+Seurat::NoLegend()+
  facet_wrap(~group, nrow=1, scales="free")

ppp=names(mn2.xdata_trajec_UDT.Int$disp)
ppp=ppp[grepl("^heat_",ppp)]
ppp=ppp[2:13]
tmpdf=data.frame()
for(x in ppp){
  qqq=mn2.xdata_trajec_UDT.Int$disp[[x]]@matrix %>% as.data.frame() %>% tidyr::pivot_longer(cols = 1:8)
  qqq$organ=gsub("heat_","",x)
  qqq$stat=sort(rep(1:5, 8))
  tmpdf=rbind(tmpdf,qqq)
}
tmpdf$name=factor(tmpdf$name, levels = sort(unique(tmpdf$name)))
tmpdf$organ=factor(tmpdf$organ, levels = unique(tmpdf$organ))
tmpdf$stat=factor(tmpdf$stat, levels = 1:5)
ggplot(tmpdf)+geom_bar(aes(organ, value, fill=stat), stat = "identity", position = position_stack())+facet_wrap(~name, nrow = 1) +xlab("")+ylab("")+theme_cowplot()+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))  ##

plot13=roveRun(xdata, Group="Location_short", CT="MPsubtype00", outdir=paste0(dir_img,"MPsub1114."), flip = FALSE, mywidth=4.5, nCT=8)

Idents(xdata)="MPsubtype00"
xdata$Cellname=rownames(xdata@meta.data)
xdata$Location_short=droplevels(xdata$Location_short)
xdata_sls_res1113=slingRun(xdata, samplesize=0.3, CTname="MPsubtype00", cellid="Cellname",
                        start.type=c("c1_CD14+.cMo"),
                        mycolor89=NULL, savesce=paste0(dir_img,"/MPsubtype00_slingshot.1113.rds"))

VlnPlot(xdata, features = c("Mono_score.01","MoDM_score.01","FeDM_score.01"), group.by = "MPsubtype00", pt.size = 0, stack = T, cols = scPalette2(30), fill.by = "ident")+geom_boxplot(width=0.3, fill="white", outlier.alpha = 0)+NoLegend()+xlab("")+ylab("")
