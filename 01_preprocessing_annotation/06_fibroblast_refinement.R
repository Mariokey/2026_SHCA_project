# ============================================================================
# Fibroblast subtype refinement, marker re-annotation and trajectory output
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L6251-6792)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
OutputDir="./AnaRes20240802/"
dir.create(OutputDir, recursive = T)

### inobj (inobj=subdataList$strom_fb)
#NRG1+:泌尿生殖系统                       FB01
#SFPR2、PLA2G2A、NT5E、NOX4:循环系统     FB02/05/09;MSC01
#PDE10A:皮肤、泌尿生殖系统              FB04
#KLF4、ADAM28:消化系统            FB08/10
#MSC02:只在前中肠出现的间质干细胞

inobj=fastRefine(inobj, dims = 1:50, resol = 0.1, nnei = 100, mindist =0.001, spread = 2)
inobj <- FindClusters(inobj, resolution = 0.01, algorithm = 2, verbose = FALSE)
DimPlot(inobj, group.by = "FB_SubClass00", reduction = "umap",label = T)
DimPlot(inobj, group.by = "seurat_clusters", reduction = "umap",label = T)

Idents(inobj)="FB_SubClass00"
xxx=FindAllMarkers(inobj, group.by = "FB_SubClass00", test.use = "t", max.cells.per.ident = 2000,
                   only.pos = TRUE, min.pct = 0.2)
AllDEG2CorPlot(inobj, xxx, identX = "FB_SubClass00",topn = 10, log_weight = 0)

Ref_UniFB_df=read.xlsx("./RefPaper_FB2021.xlsx")  ###C8-top20
Ref_UniFB=head(Ref_UniFB_df[Ref_UniFB_df$cluster==8,]$gene, n=20)

inobj=AddModuleScore(inobj, features = list(Ref_UniFB),
                     ctrl = 100, name = "UniFB0")
xxx=VlnPlot(inobj, features = "UniFB01", group.by = "FB_SubClass00", pt.size = 0)+NoLegend()
xxx_tmp=xxx[[1]]$data
xxx_tmp=xxx_tmp %>% group_by(ident) %>% summarise_at("UniFB01",median)
xxx_tmp=xxx_tmp[order(xxx_tmp$UniFB01, decreasing = TRUE),]
xxx[[1]]$data$ident=factor(xxx[[1]]$data$ident, levels = xxx_tmp$ident)
xxx+xlab("")+ggtitle("Univeral FB signature scores [Ref.]")

UniFBgenes=c("PI16","CD34")
xxx=VlnPlot(inobj, features = "PI16", group.by = "FB_SubClass00", pt.size = 0)+NoLegend()
xxx_tmp=xxx[[1]]$data
xxx1_tmp=xxx_tmp %>% group_by(ident) %>% summarise_at("PI16",mean)
xxx1_tmp=xxx1_tmp[order(xxx1_tmp$PI16, decreasing = TRUE),]
xxx=VlnPlot(inobj, features = "CD34", group.by = "FB_SubClass00", pt.size = 0)+NoLegend()
xxx_tmp=xxx[[1]]$data
xxx2_tmp=xxx_tmp %>% group_by(ident) %>% summarise_at("CD34",mean)
xxx2_tmp=xxx2_tmp[order(xxx2_tmp$CD34, decreasing = TRUE),]
xxx12=left_join(xxx1_tmp,xxx2_tmp,"ident")

inobj$cellname=rownames(inobj@meta.data)
top10 <- topDEGrun(FB_SubClass00_DEGs,topn=5)
subobj <- subset(inobj, cellname %in%
                   subsetID2(inobj@meta.data[,c('FB_SubClass00','cellname')],expected.cell = 100)$Cellname)
subobj <- ScaleData(subobj, features = rownames(subobj))
Idents(subobj) <- factor(Idents(subobj), levels = levels(inobj$FB_SubClass00))
top10$cluster <- factor(top10$cluster, levels = levels(inobj$FB_SubClass00))
top10 <- top10[order(top10$cluster),]
p1 <- DoHeatmap(subobj, features = top10$gene, size = 3, group.by = 'FB_SubClass00',angle = 30,group.bar.height = 0.05)
p1[[1]]$layers[[2]]$show.legend=F
p=p1[[1]]+scale_fill_gradient2(low = "#141464", mid = "white", high = "#78050F", midpoint = 0, limits=c(-2.5,2.5))

inobj$Organ=inobj@meta.data[,"Location_short"]
for(i in 1:length(refine_list)){
  for(j in 1:length(refine_list[[i]])){
    inobj=modifyAno(inobj, anoCol = "Organ", oldID = paste0("^",refine_list[[i]][j],"$"), newID = names(refine_list[i]))
  }
}
inobj$Organ=factor(inobj$Organ, levels = rev(names(refine_list)))
inobj$Organ=droplevels(inobj$Organ)

inobj$UMAP_1=inobj@reductions$umap@cell.embeddings[,1]
inobj$UMAP_2=inobj@reductions$umap@cell.embeddings[,2]
gpColor=hue_pal()(length(levels(inobj$Organ)))
names(gpColor)=levels(inobj$Organ)
plotlist=list()
for(x in names(gpColor)){
  ppp=ggplot()+stat_density_2d(data=inobj@meta.data %>% sample_n(size=10000),
                               aes(UMAP_1, UMAP_2, color=as.factor("..level..")), breaks=c(0.0005,100), h = c(3, 3), n=100)+
    scale_color_manual(values = c("grey",rep("NA",100),"NA"))+
    stat_density_2d(data=inobj@meta.data[inobj@meta.data$Organ %in% x, ],
                    aes_string(x = "UMAP_1", y = "UMAP_2",  alpha = "..level..", fill = "Organ"),
                    linewidth = 0, geom = "polygon",  n = 200, h = c(0.8, 0.8)) + theme_cowplot() +
    scale_fill_manual(values = gpColor)+
    scale_alpha_continuous(limits=c(0,0.3)) + NoLegend()
  plotlist[[x]]=ppp
}
plot_grid(plotlist = plotlist, ncol = 5, labels = names(plotlist))

inobj$FB_SubClass01=inobj$FB_SubClass00
rename_list=list()
rename_list[["FB01_PI16"]]=c("FB05_PLA2G2A","MSC02_CLEC3B") #,"FB09_NOX4","FB02_SFPR2","FB08_KLF4","FB13_CTRB1")
rename_list[["FB02_SFPR2"]]=c("FB02_SFPR2","FB08_KLF4", "FB09_NOX4", "FB13_CTRB1")
rename_list[["FB03_PDE10A"]]=c("FB04_PDE10A")
rename_list[["FB04_CFD"]]=c("FB03_CFD")
rename_list[["FB05_TNFAIP6"]]=c("FB06_TNFAIP6")
rename_list[["FB06_STAR"]]=c("FB07_STAR")
rename_list[["FB07_NRG1"]]=c("FB01_NRG1")
rename_list[["FB08_ADAM28"]]=c("FB10_ADAM28")
rename_list[["FB09_MMP11"]]=c("FB11_MMP11")
rename_list[["FB10_G0S2"]]=c("FB12_G0S2")
rename_list[["MSC01_NT5E"]]=c("MSC01_NT5E")
rename_list[["MyoFB01_ELN"]]=c("MyoFB01_ELN")
rename_list[["lipoFB01_APOE"]]=c("lipoFB01_APOE")
#"FB01_NRG1"     "FB02_SFPR2"    "FB03_CFD"      "FB04_PDE10A"   "FB05_PLA2G2A"
#"FB06_TNFAIP6"  "FB07_STAR"     "FB08_KLF4"     "FB09_NOX4"     "FB10_ADAM28"
#"FB11_MMP11"    "FB12_G0S2"     "FB13_CTRB1"    "MSC01_NT5E"    "MSC02_CLEC3B"
#"MyoFB01_ELN"   "lipoFB01_APOE"
for(i in 1:length(rename_list)){
  for(j in 1:length(rename_list[[i]])){
    inobj=modifyAno(inobj, anoCol = "FB_SubClass01", oldID = rename_list[[i]][j], newID = names(rename_list[i]))
  }
}

Idents(inobj)="FB_SubClass01"
FB_SubClass01_DEGs=FindAllMarkers(inobj, group.by = "FB_SubClass01", test.use = "t", max.cells.per.ident = 2000,
                                  only.pos = TRUE, min.pct = 0.2)
FB_newCTenrichFun=enrichFUNgo(ud_table = FB_SubClass01_DEGs, sp="hsa",intermode = FALSE, enriches = c("GO"))
write.xlsx(FB_newCTenrichFun$dflist$GO, file="./AnaRes20240802/FB_newCTenrichFun.GO.xlsx")
write.xlsx(FB_SubClass01_DEGs, file="./AnaRes20240802/FB_newCT_DEGs.xlsx")

inobj$cellname=rownames(inobj@meta.data)
top10 <- topDEGrun(FB_SubClass01_DEGs,topn=6)
subobj <- subset(inobj, cellname %in%
                   subsetID2(inobj@meta.data[,c('FB_SubClass01','cellname')],expected.cell = 200)$Cellname)
subobj <- ScaleData(subobj, features = rownames(subobj))
Idents(subobj) <- factor(Idents(subobj), levels = levels(inobj$FB_SubClass01))
top10$cluster <- factor(top10$cluster, levels = levels(inobj$FB_SubClass01))
top10 <- top10[order(top10$cluster),]
p1 <- DoHeatmap(subobj, features = top10$gene, size = 3, group.by = 'FB_SubClass01',angle = 30,group.bar.height = 0.05)
p1[[1]]$layers[[2]]$show.legend=F
p=p1[[1]]+scale_fill_gradient2(low = "#141464", mid = "white", high = "#78050F", midpoint = 0, limits=c(-2.5,2.5))

roveRun(inobj, Group="Organ", CT="FB_SubClass01", outdir="./AnaRes20240802/", flip = FALSE, mywidth=3.5, nCT=20)

Idents(inobj)="FB_SubClass01"
inobj$Cellname=rownames(inobj@meta.data)
inobj$Organ=droplevels(inobj$Organ)
FB_sls_res0802=slingRun(inobj, samplesize=0.2, CTname="FB_SubClass01", cellid="Cellname",
                           start.type=c("MSC01_NT5E"),
                           mycolor89=NULL, savesce=paste0("./AnaRes20240802/","/FB_inobj.slingshot.0802.rds"))
sce=readRDS(paste0("./AnaRes20240802/","/FB_inobj.slingshot.0802.rds"))
FB0802_tradeResList=tradeS(sce, samplesize=1, SvE=FALSE,  Asso=TRUE, outdirprefix=paste0("./AnaRes20240802","/FB_inobj0.1.tradeseq.0802"))
FB_sls_res0802$UMAP_allCurve+scale_size_manual(values = c(2,1.5,1,0.75,0.5,0.25))
slingshot::slingLineages(sce)

plotdf=FB_sls_res0802$plot_data
i=1
plotdf$Pseudotime=as.numeric(rescale(plotdf[,paste0("slingPseudotime_",i)], to = c(0,100)))
plotdf$Lineage=rep(paste0("Lineage", i), dim(plotdf)[1])
p2=ggplot(plotdf)+geom_density(aes(Pseudotime, color=Organ), adjust=2)+theme_cowplot()+theme(legend.position = "bottom")
p2=p2+scale_y_continuous(trans = "log1p")
zzz=ggplot(plotdf)+geom_point(aes_string("Lineage","Pseudotime", color=CTname), position = position_jitter())+
  scale_color_manual(values = mycolor89)+guides(color = guide_legend(title = "Celltype"))+
  theme_minimal_grid()+coord_flip()+theme(legend.position = "top", axis.title = element_blank(), axis.text = element_blank())+NoLegend()
mmm0=plot_grid(zzz,p2, nrow = 2,rel_heights = c(0.2,1), align = "v")
sls_res_list[[paste0("densityPoint_Curve", i)]]=mmm0

figfix="_FB"
tmpTrade=FB0802_tradeResList
GSEdir=paste0("./AnaRes20240802/","GSE_pathway.Early5Late5",figfix,"/")
dir.create(GSEdir, recursive = T)
LLLnames=names(tmpTrade)[grepl("_smoothExp$",names(tmpTrade))]
LLLnames=gsub("_smoothExp$","",LLLnames)
for(Lidx in 1:length(LLLnames)){
  LgseRes=paste0(LLLnames[Lidx],"_GSEres")
  LgseSmooth=paste0(LLLnames[Lidx],"_smoothExp")
  original_gene_list=rownames(tmpTrade[[LgseSmooth]])
  gseDF=tmpTrade[[LgseRes]]$gse_fdf
  gseOBJ=tmpTrade[[LgseRes]]$gse_obj
write.xlsx(tmpTrade[[LgseSmooth]][,1:3] %>% as.data.frame(), rowNames=TRUE, file = paste0(GSEdir,"/L",Lidx,".smoothHeatmap_col123.xlsx"))
write.xlsx(gseDF, file = paste0(GSEdir,"/L",Lidx,".GSE_pathway.EarlyLate.forDisplay.xlsx"))

toptail5=5
gseDF=gseDF[c(1:toptail5,(dim(gseDF)[1]-toptail5+1):dim(gseDF)[1]),]

for(rrr in 1:dim(gseDF)[1]){
  gseDF_row=gseDF[rrr,]
  name=gseDF_row$Description #"keratinization"
  name_go=gseDF_row$ID
  top_or_down=if_else(gseDF_row$NES>0, "head", "tail")
  show_stat=paste0("NES:\n",round(gseDF_row$NES, digit=2),"\n","p.adjust:\n",round(gseDF_row$p.adjust, digit=4))
  if(gseDF_row$pvalue<=0.05){
    JAKplot=gseaplot(gseOBJ, by = "all", title = name, geneSetID = name_go)
    runFig(paste0(GSEdir,"/L",Lidx,".NES_",gseDF_row$NES,".",gsub("\\/", ".", gsub(" ","_",name)),""), JAKplot, 6, 4.5)
  }
}
}

#split refine of FB slingshot 0811
OutputDir="./AnaRes20240812/"
dir.create(OutputDir, recursive = T)

FB_sls_res0812=slingRun(inobj,
                        samplesize=1, CTname="FB_SubClass01", cellid="Cellname",
                        start.type=c("MSC01_NT5E"),
                        mycolor89=NULL, savesce=paste0(OutputDir,"/FB_inobj.slingshot.0812.rds"))
sce=readRDS(paste0(OutputDir,"/FB_inobj.slingshot.0812.rds"))

inobj$FB_SubClass02=inobj$FB_SubClass01
inobj=modifyAno(inobj, anoCol = "FB_SubClass02", oldID = "lipoFB01_APOE", newID = "FB11_APOE")
inobj=modifyAno(inobj, anoCol = "FB_SubClass02", oldID = "FB02_SFPR2", newID = "FB02_SFRP2")
Idents(inobj)="FB_SubClass02"

calcu_list0=list()
calcu_list0[["ct_vector"]]=c("MSC01","FB01","FB02","FB05","FB11","FB07","FB08")
calcu_list0[["organ_Int"]]=c("Intestine")
calcu_list1=list()
calcu_list1[["ct_vector"]]=c("MSC01","FB01","FB02","FB05","FB11","FB07")
calcu_list1[["organ_BS"]]=c("BS")
calcu_list1[["organ_pancr"]]=c("Pancreas") #all are FB07_NRG1 (40 cells only)
calcu_list1[["organ_UDT"]]=c("UDT")
calcu_list1[["organ_UB"]]=c("UB")
calcu_list1[["organ_US"]]=c("Kidney") # run wrong! 66 cells only
calcu_list1[["organ_VDFU"]]=c("VD","FU")
calcu_list1[["organ_TesOva"]]=c("Testis","Ovary")
calcu_list2=list()
calcu_list2[["ct_vector"]]=c("MSC01","FB01","FB02","FB05","FB11","FB03")
calcu_list2[["organ_VDFU"]]=c("VD","FU")
calcu_list2[["organ_TesOva"]]=c("Testis","Ovary")
calcu_list2[["organ_Skin"]]=c("Skin")
calcu_list3=list()
calcu_list3[["ct_vector"]]=c("MSC01","FB01","FB02","FB05","FB11","FB06")
calcu_list3[["organ_liver"]]=c("Liver") #34 cells
calcu_list3[["organ_pancr"]]=c("Pancreas") #9 cells
calcu_list3[["organ_US"]]=c("Kidney") #41 cells
calcu_list3[["organ_VDFU"]]=c("VD","FU")
calcu_list3[["organ_TesOva"]]=c("Testis","Ovary")

for(idx in 1:3){
  idx=3
  mark=paste0("L",idx)
  tmplist=eval(parse(text=paste0("calcu_list",idx)))

sub_inobj=subset(inobj, FB_SubClass02 %in%
                   lapply(tmplist[[1]], FUN = function(i)levels(inobj$FB_SubClass02)[grep(paste0("^",i), levels(inobj$FB_SubClass02))]) %>% unlist())
sub_inobj$FB_SubClass02=droplevels(sub_inobj$FB_SubClass02)
x=4
basemark=paste0(mark,"_",names(tmplist)[x],"_")
FBres_traj[[paste0(basemark,"seudata")]]=subset(sub_inobj, Organ %in% tmplist[[x]])
FBres_traj[[paste0(basemark,"seudata")]]$Cellname=rownames(FBres_traj[[paste0(basemark,"seudata")]]@meta.data)
FBres_traj[[paste0(basemark,"seudata")]]$Organ=droplevels(FBres_traj[[paste0(basemark,"seudata")]]$Organ)
FBres_traj[[paste0(basemark,"seudata")]]=fastRefine(FBres_traj[[paste0(basemark,"seudata")]],
                                                    dims = 1:50, resol = 0.1,
                                                    nnei = 20,
                                                    mindist =0.35, spread = 0.5)
FBres_traj[[paste0(basemark,"slingshot")]]=slingRun(FBres_traj[[paste0(basemark,"seudata")]],
                                                    samplesize=1,
                                                    start.type=c("FB01_PI16"), #c("MSC01_NT5E"),
                      CTname="FB_SubClass02", cellid="Cellname",
                      mycolor89=NULL, savesce=paste0(OutputDir,"/",paste0(basemark,"slingshot"),".rds"))
sce=readRDS(paste0(OutputDir,"/",paste0(basemark,"slingshot"),".rds"))
FBres_traj[[paste0(basemark,"tradeS")]]=tradeS(sce,
                                               samplesize=1,
                                               SvE=FALSE,  Asso=TRUE,
                                outdirprefix=paste0(OutputDir,"/",paste0(basemark,"trades")))

FB_sls_res=FBres_traj[[paste0(basemark,"slingshot")]]
max_lin=max(FB_sls_res$plot_curve_data$Lineage)
xxx=FB_sls_res$UMAP_allCurve+scale_size_manual(values = rev(1:max_lin %>% rescale(to=c(0.25,2))))
runFig(paste0(OutputDir,"/",paste0(basemark,"sling_UMAP")), xxx, 6, 4.5)
for(xi in 1:max_lin){
xxx=FB_sls_res[[paste0("densityPoint_Curve",xi)]]
runFig(paste0(OutputDir,"/",paste0(basemark,"sling_density_",xi)), xxx, 6, 4.5)
}

tmpTrade=FBres_traj[[paste0(basemark,"tradeS")]]
GSEdir=paste0(OutputDir,"/",paste0(basemark,"tradeS"),"_GSE_pathway.Early5Late5/")
dir.create(GSEdir, recursive = T)
LLLnames=names(tmpTrade)[grepl("_smoothExp$",names(tmpTrade))]
LLLnames=gsub("_smoothExp$","",LLLnames)
for(Lidx in 1:length(LLLnames)){
  LgseRes=paste0(LLLnames[Lidx],"_GSEres")
  LgseSmooth=paste0(LLLnames[Lidx],"_smoothExp")
  original_gene_list=rownames(tmpTrade[[LgseSmooth]])
  gseDF=tmpTrade[[LgseRes]]$gse_fdf
  gseOBJ=tmpTrade[[LgseRes]]$gse_obj
write.xlsx(tmpTrade[[LgseSmooth]][,1:3] %>% as.data.frame(), rowNames=TRUE,
           file = paste0(GSEdir,"/L",Lidx,".smoothHeatmap_col123.xlsx"))
write.xlsx(gseDF, file = paste0(GSEdir,"/L",Lidx,".GSE_pathway.EarlyLate.forDisplay.xlsx"))

toptail5=5
gseDF=gseDF[c(1:toptail5,(dim(gseDF)[1]-toptail5+1):dim(gseDF)[1]),]

for(rrr in 1:dim(gseDF)[1]){
  gseDF_row=gseDF[rrr,]
  name=gseDF_row$Description #"keratinization"
  name_go=gseDF_row$ID
  top_or_down=if_else(gseDF_row$NES>0, "head", "tail")
  show_stat=paste0("NES:\n",round(gseDF_row$NES, digit=2),"\n","p.adjust:\n",round(gseDF_row$p.adjust, digit=4))
  if(gseDF_row$pvalue<=0.05){
    JAKplot=gseaplot(gseOBJ, by = "all", title = name, geneSetID = name_go)
    runFig(paste0(GSEdir,"/L",Lidx,".NES_",gseDF_row$NES,".",gsub("\\/", ".", gsub(" ","_",name)),""), JAKplot, 6, 4.5)
  }
}
}

}

### 0819 add refineUMAP/DEGs/Function for Organs (FB07_NRG need further confirm?)
resDir="./AnaRes20240819/"
dir.create(resDir, recursive = TRUE)

tslist=list()
tslist[["BS"]]=c("BS")
tslist[["UDT"]]=c("UDT")
tslist[["Intestine"]]=c("Intestine")
tslist[["LungAirway"]]=c("Lung","Airway")
tslist[["UB"]]=c("UB")
tslist[["Artery"]]=c("Artery")
tslist[["Skin"]]=c("Skin")
tslist[["VDFU"]]=c("VD","FU")
tslist[["TestisOvary"]]=c("Testis","Ovary")

FBct02_degfun_list=list()
for(x in 1:length(tslist)){
  #x=9
  xname=names(tslist)[x]
  fct=tslist[[x]]
  tmpobj=subset(inobj, Organ %in% fct)
  xtable=table(tmpobj$FB_SubClass02)
  xtable=xtable[xtable>=30]
  tmpobj=subset(tmpobj, FB_SubClass02 %in% names(xtable))
  tmpobj$FB_SubClass02=droplevels(tmpobj$FB_SubClass02)
  tmpobj$Organ=droplevels(tmpobj$Organ)
  tmpobj$Location_short=droplevels(tmpobj$Location_short)

  tmpobj=fastRefine(tmpobj, dims = 1:50, resol = 0.1, nnei = 60, mindist =0.25, spread = 0.25)
  Idents(tmpobj)="FB_SubClass02"
  xxx=DimPlot_idx(tmpobj, group.by = "FB_SubClass02",
                  cols=scPalette1(length(levels(tmpobj$FB_SubClass02))),
                  reduction = "umap", repel = T)
  runFig(paste0(resDir,"/",xname,"_UMAP"), xxx, 6.5, 4.5)
  if(length(names(table(tmpobj$Location_short)))>1){
  yyy=fastpHeat(table(tmpobj$FB_SubClass02, tmpobj$Location_short), dispNum = TRUE, rowcol = "c")
  runFig(paste0(resDir,"/",xname,"_heat_ctPctNum"), yyy, 1.5+1.5*(length(names(table(tmpobj$Location_short)))-1), 4.5)
  }
  FBct02_degfun_list[[paste0(xname,"_seudata")]]=tmpobj

  FBct02_degfun_list[[xname]]=DEG2FUN(tmpobj, myiden="FB_SubClass02", groupname="Organ",
         pval_threshold = 0.01, logfc_threshold = 0.2, min_pct = 0.1, DEGenefilter=FALSE,
         DEGct=TRUE, DEGgroup=FALSE,
         maxCellperIdents=1000,
         FunctionGO=TRUE, FunctionKEGG=TRUE, FunctionGSE=TRUE,
         OutDir = resDir, suffix = paste0(xname,"_FBct02"))
}

fastpHeat=function(zzz, dispNum=FALSE, rowcol="r"){
  if(dispNum){
    dispNum=zzz #round(zzz/1000,digits = 2)
  }
  if(rowcol=="r"){
    zzzr=rowPercents(zzz)[,1:dim(zzz)[2]]
    pheatmap(zzzr, name = "pct.",
             display_numbers = dispNum,
             cluster_rows = F, cluster_cols = F
             )
  }
  if(rowcol=="c"){
    zzzc=colPercents(zzz)[1:dim(zzz)[1],]
    pheatmap(zzzc, name = "pct.",
             display_numbers = dispNum,
             cluster_rows = F, cluster_cols = F
             )
  }
}

FBinobj=subdataList$strom_fb
fbfiles=list.files("./AnaRes20240819/")
fbfiles=fbfiles[grepl("DEGfilter",fbfiles)]
outDegList=list()
for(x in fbfiles){
  tmpDF=read.xlsx(paste0("./AnaRes20240819/",x))
  tmpDF$cluster=tmpDF$celltype
outDegList[[gsub("_.+","",x)]]=tmpDF
}
#p_val avg_log2FC pct.1 pct.2    p_val_adj   cluster     gene rows
#HLA-DQA1 8.426078e-63   2.271761 1.000 0.320 2.866215e-58 DC01_CD1C HLA-DQA1    2

### 20251130 add final output
outdir="./AnaRes20240819/"
p=scPioneer::DimPlot_idx(subdataList$strom_fb, group.by = "FB_SubClass02",
            cols=scPioneer::scPalette1(length(levels(subdataList$strom_fb$FB_SubClass02))),
            reduction = "umap", repel = T)+NoAxes()
ggsave(paste0(outdir,'SX_FB_a.',"UMAP",'.pdf'),p,width = 5.5, height = 3.6)

p=roveRun(subdataList$strom_fb, Group="Organ", CT="FB_SubClass02", outdir="./AnaRes20240819/", flip = FALSE, mywidth=5, nCT=14)

p=scPioneer::FeaturePlot2(subdataList$strom_fb, features = c("CD34","PI16","COL15A1","CXCL12"), ncols = 2)
ggsave(paste0(outdir,'SX_FB_c.',"4gene",'.pdf'),p,width = 5.5, height = 4.5)

library(CytoTRACE)
xdata_results = CytoTRACE(as.matrix(inobj@assays$RNA$counts), ncores = 4, subsamplesize = 1000)
ptV=inobj$FB_SubClass02 %>% as.character()
names(ptV)=names(inobj$FB_SubClass02)
### NOTE DO NOT load scPioneer!!! ->boxplot will be replaced by this package!!!
plotCytoTRACE(xdata_results, phenotype = ptV,  outputDir = outdir,
              emb = inobj@reductions$umap@cell.embeddings %>% as.data.frame()) #, gene = "Kit")
inobj$CytoTRACE=xdata_results$CytoTRACE
p=scPioneer::FeaturePlot2(inobj, features = c("CytoTRACE"))
ggsave(paste0(outdir,'SX_FB_d.',"CytoTRACE",'.pdf'),p,width = 3.5, height = 3)

seleFB_terms=c("axon guidance", "neuron projection guidance", "axonogenesis", "synapse organization", "forebrain development", "modulation of chemical synaptic transmission", "morphogenesis of a branching structure", "muscle tissue development", "epithelial cell proliferation", "axonogenesis", "branching morphogenesis of an epithelial tube", "mammary gland development", "response to oxidative stress", "response to reactive oxygen species", "cellular response to oxidative stress", "cellular response to reactive oxygen species", "cellular oxidant detoxification", "response to hydrogen peroxide", "alpha-beta T cell differentiation", "alpha-beta T cell activation", "CD4-positive, alpha-beta T cell differentiation", "T-helper 17 cell differentiation", "T-helper cell differentiation", "CD4-positive, alpha-beta T cell differentiation involved in immune response", "organelle disassembly", "sequestering of metal ion", "regulation of fibroblast proliferation", "detoxification", "regulation of transmembrane receptor protein serine", "threonine kinase signaling pathway", "extracellular matrix organization", "extracellular structure organization", "external encapsulating structure organization", "myoblast differentiation", "insulin-like growth factor receptor signaling pathway", "positive regulation of insulin-like growth factor receptor signaling pathway", "non-canonical Wnt signaling pathway", "positive regulation of MAPK cascade", "regulation of ERK1 and ERK2 cascade", "cell-substrate adhesion", "negative chemotaxis", "outflow tract morphogenesis", "regulation of cell-substrate adhesion", "extracellular matrix organization", "sensory organ morphogenesis")
allFB_terms=read.xlsx("./AnaRes20240913/FB_8Module.FuncEnrich.xlsx")
allFB_terms_select=allFB_terms[allFB_terms$Description %in% seleFB_terms,]
allFB_terms_select$Description=factor(allFB_terms_select$Description, levels=unique(seleFB_terms))
p=ggplot(allFB_terms_select, aes(Module, Description))+
  geom_point(aes(size=Factor,fill=-log10(pvalue)), shape=22)+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = 45)+
  xlab("")+ylab("")+
  theme_cowplot()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.y = element_blank()) #+
ggsave(paste0(outdir,'SX_FB_f.',"ModuleFun",'.pdf'),p,width = 10, height = 7)

outdir="./MP_reScore241104/"
inobj=readRDS("./MP_reScore241104/MPsubtype00_slingshot.1113.subseu.rds")
p=scPioneer::DimPlot_idx(inobj, group.by = "MPsubtype00",
                         cols=scPioneer::scPalette2(30), #scPioneer::scPalette1(length(levels(inobj$MPsubtype00))),
                         reduction = "umap", repel = T)+NoAxes()
ggsave(paste0(outdir,'SX_MP_a.',"UMAP",'.pdf'),p,width = 5, height = 3.2)

p=roveRun(inobj, Group="Organ", CT="MPsubtype00", outdir=outdir, flip = FALSE, mywidth=5, nCT=9)

p=scPioneer::FeaturePlot2(inobj, features = c("FCN1","FABP4","STARD13"), ncols = 3)
ggsave(paste0(outdir,'SX_MP_b2.',"3gene",'.pdf'),p,width = 6, height = 2.2)

p=VlnPlot(inobj, features = c("Mono_score.01","MoDM_score.01","FeDM_score.01"), group.by = "MPsubtype00", pt.size = 0, stack = T, cols = scPioneer::scPalette2(30), fill.by = "ident")+geom_boxplot(width=0.3, fill="white", outlier.alpha = 0)+NoLegend()+xlab("")+ylab("")
ggsave(paste0(outdir,'SX_MP_b1.',"3scores",'.pdf'),p,width = 7, height = 5)

mysce=readRDS("./MP_reScore241104/MPsubtype00_slingshot.1113.rds")
