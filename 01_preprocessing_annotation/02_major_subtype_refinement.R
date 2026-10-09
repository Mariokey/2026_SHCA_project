# ============================================================================
# Major class and subtype label refinement across all compartments
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L2332-2917)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### ??? alternative strategy to manipulate large dataset:
###     use seurat5 Sketch & ProjectData to analyze the data
Mobj=readRDS("./Mobj.withX_unk.v20240524.rds")
Mobj$refine_ct00=Mobj$subtypes03
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Xendo_unk.", newID = "Mast",
               newOrder = c(levels(Mobj$refine_ct00)[1:11],
                            "Mast",
                            levels(Mobj$refine_ct00)[12:51]))
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Xstrom_unk.", newID = "FB")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "MSC", newID = "FB")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Pericytes", newID = "SMC")

Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "long_ST", newID = "ST")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "round_ST", newID = "ST")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "early_SPC", newID = "SPC")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "late_SPC", newID = "SPC")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Xgerm_unk.", newID = "ST")

Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "ILC3", newID = "T")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "CD4T", newID = "T")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "CD8T", newID = "T")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "gdT", newID = "T")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "DNT", newID = "T")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "NKT", newID = "T")

Mobj$cellid_ori=Mobj$cellid
Mobj$cellid=rownames(Mobj@meta.data)
Mobj=modifyAno(Mobj, anoCol = "refine_ct00",
               cellid = colnames(subset(Mobj, subtypes01 %in% c("B03_MZB1"))),
               newID = "Plasma")

Mobj=subset(Mobj, refine_ct00 %!in% c("Xepi_unk.","Ximmu_unk."))

ori_vector=c('Basal_Cells', 'Early_Suprabasal_Cells', 'Mid_Suprabasal_Cells', 'Early_Suprabasal_Cells', 'Prezymogenic_Cells', 'Chief_Cells', 'Pit_Cells', 'Metallothionein_Cells', 'Enterocytes', 'Goblet_Cells', 'BEST4_Epithelial_Cells', 'Colonocytes', 'Club_Cells', 'Multiciliated_Cells', 'Airway_Suprabasal_Cells', 'Cortical_Collecting_Duct_Principal_Cells', 'Papillary_Epithelial_Cells', 'Melanocytes', 'Spinous_Cells', 'Vellus_Hair_Follicle', 'Hepatocytes', 'Endocrine_Cells', 'Tuft_Cells', 'DS_Epithelial_Cells_1', 'Mesothelial_Cells', 'RPS_CS_Epithelial_Cells_1', 'RPS_CS_Epithelial_Cells_2', 'MIX_Epithelial_Cells_1', 'RPS_Epithelial_Cells_2', 'Pancreas_Cells', 'Inner_Hair_Follicle')#,  'Glia_Cells', 'Neurons',  'Plasma_Cells', 'Immune_Cells')
new_vector=c('Basal', 'Superbasal', 'Superbasal', 'Superbasal', 'Gastric', 'Gastric', 'Gastric', 'Gastric', 'Intestinal', 'Intestinal', 'Intestinal', 'Intestinal', 'Tracheal', 'Tracheal', 'Tracheal', 'Collecting_ductal', 'Collecting_ductal', 'Skin', 'Skin', 'Skin', 'Hepatocyte', 'Endocrine', 'Tuft', 'BEC', 'Mesothelial', 'Epi_undef1', 'Epi_undef2', 'Epi_undef3', 'RPS_Epi_undef4', 'Pancreas', 'Skin')#, 'Glia', 'Neuron',  'Immune', 'Immune')
for(i in 1:length(ori_vector)){
  print(ori_vector[i])
Mobj=modifyAno(Mobj, anoCol = "refine_ct00",
               cellid = colnames(subset(Mobj, subtypes01 %in% ori_vector[i])), #oldID = ori_vector[i],
               newID = new_vector[i])
}
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Granulosa.Epi", newID = "Granulosa")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Collecting-ductal", newID = "Collecting_ductal")
Mobj=modifyAno(Mobj, anoCol = "refine_ct00", oldID = "Suprabasal", newID = "Superbasal")

Mobj=subset(Mobj, refine_ct00 %!in% c("Xepi_unk.","Ximmu_unk."))
CT33_order=c("B","Plasma",
             "T","NK",
             "Mono","Mph","Neutro","Mast",
             "Basal","Superbasal","Tuft","Hepatocyte","BEC","Pancreas",
             "Endocrine","Collecting_ductal","Granulosa","Mesothelial",
             "Gastric","Intestinal","Tracheal","Skin",
             "Epi_undef1","Epi_undef2","Epi_undef3","RPS_Epi_undef4",
              "EC","FB","SMC","SPG","SPC","ST","Prolif." )
Mobj$refine_ct00=factor(Mobj$refine_ct00, levels = CT33_order)
DimPlot(Mobj, group.by = "refine_ct00", reduction = "umap", label = T, repel = T)+NoLegend()

iCT=levels(Mobj$refine_ct00)
nCT=length(iCT)
mycolortmp=sample(sciColor_filter$hex_color, size = nCT)
show_col(mycolortmp)
names(mycolortmp)=iCT
DimPlot_idx(Mobj, group.by = 'refine_ct00', cols = as.vector(mycolortmp), repel = T, label.idx.size = 5, raster = T)+xlab("UMAP")+ylab("UMAP")+myfont

xxx=FindMarkers(Mobj, ident.1 = c("Xgerm_unk."), group.by = "refine_ct00", test.use = "t",
                max.cells.per.ident = 3000, only.pos = T, logfc.threshold = 0.5, min.pct = 0.3)
head(xxx, n=30)

DimPlot(Mobj, cells.highlight = colnames(subset(Mobj, refine_ct00 %in% c("granulosa"))),
        reduction = "umap", label = T, repel = T)+NoLegend()

source(file = "DATA_DIR")
CellularModules(Mobj,annotation = "refine_ct00", Sample = "sample",Disease = "Location_short",Split_num = 8,OutDir = "./CMpredict/", )

source(file = "DATA_DIR")
Mobj_tmp=subset(Mobj, refine_ct00 %in% c("B","Plasma","T","NK",
                                         "Mono","Mph","Neutro","Mast",
                                         "EC","FB","SMC"))
table(Mobj_tmp$subtypes011)
Mobj_tmp$subtypes011=droplevels(Mobj_tmp$subtypes011)
Mobj_tmp$sample=droplevels(Mobj_tmp$sample)
Mobj_tmp$Location_short=droplevels(Mobj_tmp$Location_short)
tmpmat=table(Mobj_tmp$Location_short, Mobj_tmp$sample)
tmpmat[tmpmat<=5]=0
tmpmat=tmpmat[names(rowSums2(tmpmat)>50),]
tmpmat=tmpmat[,names(colSums2(tmpmat)>50)]
CellularModules(Mobj_tmp,annotation = "subtypes011", toptext_size1 = 30, righttext_size1 = 30, Sample = "sample",Disease = "Location_short",Split_num = 5,OutDir = "./CMpredict_IMStroma/")

source(file = "DATA_DIR")
CellularModules(Mobj_keep,annotation = "refine_ct03", Sample = "sample",Disease = "Location_short",Split_num = 8,OutDir = "./CMpredict_ct2502/")
Mobj_tmp=subset(Mobj_keep, refine_ct03 %!in% c("Prolif.","ST","SPC","SPG"))
table(Mobj_tmp$subtypes015)
Mobj_tmp$subtypes015=droplevels(Mobj_tmp$subtypes015)
Mobj_tmp$sample=droplevels(Mobj_tmp$sample)
Mobj_tmp$Location_short=droplevels(Mobj_tmp$Location_short)
tmpmat=table(Mobj_tmp$Location_short, Mobj_tmp$sample)
tmpmat[tmpmat<=5]=0
tmpmat=tmpmat[names(rowSums(tmpmat)>50),]
tmpmat=tmpmat[,names(colSums(tmpmat)>50)]
CellularModules(Mobj_tmp,annotation = "subtypes015", toptext_size1 = 30, righttext_size1 = 30, Sample = "sample",Disease = "Location_short",Split_num = 8,OutDir = "./CMpredict_sub2502/")

#xxx=read.csv("./CMpredict_IMStroma/2.CM_Composition.composition.csv")
xxx=read.csv("./CMpredict_sub2502/2.CM_Composition.composition.csv")
xxx_long=melt(xxx)
colnames(xxx_long)=c("Module","Sample","Value")
xxx_long$System=gsub("_.+","",xxx_long$Sample)
xxx_long$Location=gsub(".+_","",xxx_long$Sample)
xxx_long$Location=gsub("\\d+","",xxx_long$Location)
xxx_long$Location=gsub("^TS$","Testis",xxx_long$Location)
xxx_long$Location=gsub("^LMB$","MB",xxx_long$Location)
xxx_long$Location=gsub("^RMB$","MB",xxx_long$Location)
xxx_long$Location=factor(xxx_long$Location, levels = levels(Mobj_keep$Location_short))
xxx_long$Module =factor(xxx_long$Module , levels = paste0("CM",1:8))
head(xxx_long)
ggplot()+
  geom_boxplot(data=xxx_long, aes(Location, Value), width=0.5, fill="white", outlier.alpha = 0)+xlab("")+theme_cowplot()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+facet_wrap(~Module, ncol = 1)

yyy=read.table("./rove_main_age.txt", sep = "\t", header = T)
yyy$Celltype=rownames(yyy)
yyy_long=melt(yyy,id.vars = "Celltype")
colnames(yyy_long)=c("Celltype","Location","Roe")
yyy_long$Location=factor(yyy_long$Location, levels = loca_order)
yyy_long$Celltype=factor(yyy_long$Celltype, levels = levels(Mobj$subtypes011))
ggplot()+geom_bar(data=yyy_long,aes(Location, Roe, fill=Celltype), stat = "identity")
ggplot()+geom_bar(data= yyy_long, #%>% group_by(Location) %>% top_n(3,Roe),
                  aes(Location, Roe, fill=Celltype), stat = "identity")+theme_minimal_grid()+coord_polar()+NoLegend()

DEGS=readRDS(file =paste0("./CMpredict/","/","3.CM_ExpressionRadar_DEG.rds"))
DEGS=DEGS[!grepl("^A[CPL]\\d+|-AS1$|^LINC",DEGS$gene),]
DEGS %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC*(pct.1-pct.2)) -> top5
CMgenelist <- split(top5$gene, top5$cluster)

Outgene <- unlist(CMgenelist)#[!unlist(CMgenelist) %in% rownames(OBJ_CMS@assays[[assay]])]

ring_start.degree  = 90+180/length(CMgenelist)
Rader_ring(OBJ_CMS, assay = assay, CM_metaname = CM_metaname, genelist = CMgenelist, ColorCMS = ColorCMS,
         ring_start.degree = ring_start.degree, #-360/2/5, #155-10*(Split_num - 2),
         gene_labelcex = Radar_labelcex,
         Filename = paste0(OutDir,"/","3.CM_ExpressionRadar") )

mat_2=readRDS("./CMpredict/3.CM_ExpressionRadar_tmp_mat_2.rds")
ggradar(mat_2,
        group.point.size = 5,
        gridline.mid.colour = "#4097a0", gridline.max.colour = "grey",grid.line.width = 1,group.line.width = 3,
        axis.label.size = 5, axis.line.colour = "white",
        plot.legend = F, # legend.text.size = 14, legend.position = "left",
        background.circle.colour = "white",
        background.circle.transparency = 0,
        label.gridline.max = F, label.gridline.min = F, label.gridline.mid = F,
        group.colours = pal_npg(alpha = 0.9)(10)[1:8])

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
colnames(yyy)=c("group","CT","proportion")
yyy$CT=factor(yyy$CT, levels = rev(levels(Mobj_tmp$subtypes011)))
yyy$group=factor(yyy$group, levels = levels(Mobj_tmp$System))

zzz1=yyy
zzz1$group=if_else(zzz1$group %in% c("CS"),"CS", "other")
addDF1=aggregate(proportion ~ group + CT, data = zzz1, FUN = sum, na.rm = TRUE)
addDF1=addDF1[addDF1$group %in% c("CS"),]
zzz2=yyy
zzz2$group=if_else(zzz2$group %in% c("CS","DS"),"DS", "other")
addDF2=aggregate(proportion ~ group + CT, data = zzz2, FUN = sum, na.rm = TRUE)
addDF2=addDF2[addDF2$group %in% c("DS"),]
zzz3=yyy
zzz3$group=if_else(zzz3$group %in% c("CS","DS","RPS"),"RPS", "other")
addDF3=aggregate(proportion ~ group + CT, data = zzz3, FUN = sum, na.rm = TRUE)
addDF3=addDF3[addDF3$group %in% c("RPS"),]
zzz4=yyy
zzz4$group=if_else(zzz4$group %in% c("CS","DS","RPS","RS"),"RS", "other")
addDF4=aggregate(proportion ~ group + CT, data = zzz4, FUN = sum, na.rm = TRUE)
addDF4=addDF4[addDF4$group %in% c("RS"),]
zzz5=yyy
zzz5$group=if_else(zzz5$group %in% c("CS","DS","RPS","RS","SKIN"),"SKIN", "other")
addDF5=aggregate(proportion ~ group + CT, data = zzz5, FUN = sum, na.rm = TRUE)
addDF5=addDF5[addDF5$group %in% c("SKIN"),]

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
scale_fill_manual(values = hue_pal()(6))+#+theme_minimal()+theme(axis.text = element_blank())+ #theme(axis.text = element_text(size=14, colour = "black"))+
xlab("")+ylab("")+coord_flip()#+NoLegend()
runFig(paste0(outdir,"barplot.6Sys105CT_NoMastProlif"), yyy_plot, 5.75, 13)

ncols=1
xxxplot=DimPlot_idx(Mobj_tmp, group.by = "subtypes011", idx.sep = "  ", prefix.index = NULL,
                    cols=as.vector(mycolortmp[names(mycolortmp) %!in% c("Mast","Prolif.")]))+
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
    lll[[li]][["size"]]=unit(0.675, "cm")
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
    label_col <- ifelse(hcl[, "l"] > 50, "black", "white")
    lll[[li]]$children[[1]]$children[[1]]$gp$col=as.vector(label_col)
  }
}
yyy[["grobs"]][[1]][["grobs"]]=lll
Bplot=plot_grid(NULL,
                plot_grid(NULL,yyy,nrow = 3, rel_heights = c(0.01,0.98,0.01)),
                NULL,
                nrow = 1, rel_widths = c(0.12,0.87,0.01))
ppp=Aplot+annotation_custom(grob = ggplotGrob(Bplot))
runFig(paste0(outdir, "/All105.6Sys105CT_NoMastProlif.CTlegend.test"), ppp, 4.25, 22.5)

roveRun(Mobj, Group="Location_short", CT="subtypes011", outdir="./", mywidth=20, nCT=107)

################################################################################
################ Color for current project 20240619  ##################

pick_group_color=function(color_df, num_clusters=24, col2number=list(ct=c("A","B"),num=c(2,3)), ...){

  color_res=list()
  # Number of clusters/groups you want

  # Perform k-means clustering on the RGB values
  set.seed(42) # for reproducibility
  kmeans_result <- kmeans(colors_df[, c("r", "g", "b")], centers = max(num_clusters,length(col2number[["ct"]])))

  # Add the cluster assignment to the data frame
  colors_df$cluster <- kmeans_result$cluster

  # Function to get colors belonging to each group
  get_color_group <- function(cluster_number, df, n_colors = 30) {
    group_colors <- df %>%
      filter(cluster == cluster_number) %>%
      slice_head(n = n_colors) # Take the top n colors
    return(group_colors)
  }

  for(i in 1:length(col2number[["ct"]])){
  ## Example: Get colors from group 1
    color_res[[col2number[["ct"]][i]]] <- get_color_group(i, colors_df, col2number[["num"]][i])$hex_color
  }

  return(color_res)
}

### 20240831
outdir="./AnaRes20240829/"
Mobj$MajorClass02=Mobj$MajorClass01
Mobj=modifyAno(Mobj, anoCol = "MajorClass02",
               cellid = colnames(subset(Mobj, refine_ct00 %in% "Mast")),
               newID = "Immune Cells",
               newOrder = c("Immune Cells","Epithelial Cells","Endothelial Cells","Stromal Cells","Germ Cells","Prolif."))

Mobj$Organ=Mobj@meta.data[,"Location_short"]
for(i in 1:length(refine_list)){
  for(j in 1:length(refine_list[[i]])){
    Mobj=modifyAno(Mobj, anoCol = "Organ", oldID = paste0("^",refine_list[[i]][j],"$"), newID = names(refine_list[i]))
  }
}
Mobj$Organ=factor(Mobj$Organ, levels = rev(names(refine_list)))
Mobj$Organ=droplevels(Mobj$Organ)

roveRun(Mobj, Group="Organ", CT="refine_ct00", outdir=outdir, flip = FALSE, mywidth=5, nCT=31)
roveRun(Mobj, Group="Location_short", CT="refine_ct00", outdir=outdir, flip = FALSE, mywidth=5, nCT=31)

metaf <- Mobj@meta.data[,c('Organ',"Organ","refine_ct00")]
p <- scPioneer::plotbox(metaf)
ggsave(paste0(outdir,'plotFracBox_splitby_',"Organ",'.pdf'),p,width = 12*(31/23), height = 4.5)

metaf <- Mobj@meta.data[,c("Organ","refine_ct00")]
tmpCTnums=colSums(table(Mobj$refine_ct00, Mobj$MajorClass02)>0)
major2ct_color_list=pick_group_color(sciColor_filter, num_clusters = 8,
                 col2number = list(ct=names(tmpCTnums), num=as.vector(tmpCTnums)))
CTcolor=unlist(major2ct_color_list)
names(CTcolor)=levels(Mobj$refine_ct00)
p3 <- plot_fraction(metaf, color.use = CTcolor)
ggsave(paste0(outdir, 'Barplot_',"Organ",'.pdf'), p3, width = 9*(1+0.2), height = 7)
p4 <- plot_fraction(obj@meta.data[,c('Sample',CT)], color.use = CTcolor)
ggsave(paste0(outdir, 'Barplot_Sample.pdf'), p4, width = 9*(1+addwidth), height = 7)

DimPlot_idx(Mobj, group.by = "refine_ct00", cols=as.vector(CTcolor))
DimPlot_idx(Mobj, group.by = "System", label.idx = FALSE)+NoAxes()
DimPlot_idx(Mobj, group.by = "Organ", label.idx = FALSE)+NoAxes()
DimPlot_idx(Mobj, group.by = "Location_short", label.idx = FALSE)+NoAxes()
DimPlot_idx(Mobj, group.by = "Donor", label.idx = FALSE)+NoAxes()

xxx=VlnPlot(Mobj, features = c("nCount_RNA","nFeature_RNA","pct_counts_mt"),
            group.by = "sample", pt.size = 0, ncol = 1)+geom_boxplot(width=0.25)+NoLegend()+xlab("")
xxx[[1]]=xxx[[1]]+geom_boxplot(width=0.25)
xxx[[2]]=xxx[[2]]+geom_boxplot(width=0.25)
xxx= vln_modif(xxx)
ggsave(paste0(outdir, 'vlnplot_QCbySample.pdf'), xxx, width = 24, height = 7.5)

xxx=VlnPlot(Mobj, features = c("nCount_RNA","nFeature_RNA","pct_counts_mt"),
            group.by = "Donor", pt.size = 0, ncol = 1)+geom_boxplot(width=0.25)+NoLegend()+xlab("")
xxx[[1]]=xxx[[1]]+geom_boxplot(width=0.25)
xxx[[2]]=xxx[[2]]+geom_boxplot(width=0.25)
xxx= vln_modif(xxx)
ggsave(paste0(outdir, 'vlnplot_QCbyDonor.pdf'), xxx, width = 15, height = 7.5)

seu_immu$forCombine %>% levels()
seu_epi$forCombine %>% levels()   #
seu_endo$forCombine %>% levels()  #"Xendo_unk."
seu_strom$forCombine %>% levels() #"Xstrom_unk."
seu_germ$forCombine %>% levels()  #"granulosa"  "Xgerm_unk."
#xxx=DimPlot_idx(, group.by = "CT89", prefix.index = 28:37, cols = as.vector(CT89_color[levels(dData_SgMelSch$CT89)]), raster = T, raster.dpi=c(1024*0.4,1024*0.4), repel =T)+NoLegend()

old_seu=c("seu_immu","seu_epi","seu_endo","seu_strom","seu_germ")
All107CTs=levels(Mobj$subtypes011)
All107CTs_1=All107CTs[1:31]
All107CTs_2=All107CTs[32:61]
All107CTs_3=All107CTs[62:72]
All107CTs_4=All107CTs[73:97]
All107CTs_5=All107CTs[98:106]
reSeuList=list()
for(x in 1:5){
  #x=1
  xxx=eval(parse(text = paste0("All107CTs_",x)))
  old_data=subset(Mobj, subtypes011 %in% xxx)
  old_data$subtypes011=droplevels(old_data$subtypes011)
  new_data=eval(parse(text = old_seu[x]))
  new_data$cellid=rownames(new_data@meta.data)
  new_data=subset(new_data, cellid %in% colnames(old_data))
  new_data$subtypes00=new_data$forCombine
  for(xct in levels(old_data$subtypes011)){
  new_data=modifyAno(new_data, anoCol = "subtypes00",
                 cellid = colnames(subset(old_data, subtypes011 %in% xct)),
                 newID = xct,
                 newOrder = levels(old_data$subtypes011))
  }
  reSeuList[[old_seu[x]]]=new_data

}

#Mobj$subtypes011
iCT=levels(Mobj$subtypes011)
nCT=length(iCT)
mycolortmp=sample(sciColor_filter$hex_color, size = nCT)
show_col(mycolortmp)
names(mycolortmp)=iCT

DimPlot_idx(reSeuList$seu_immu, group.by = "subtypes00", prefix.index = 1:30,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_immu$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
DimPlot_idx(reSeuList$seu_epi, group.by = "subtypes00", prefix.index = 31:60,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_epi$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
DimPlot_idx(reSeuList$seu_endo, group.by = "subtypes00", prefix.index = 61:71,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_endo$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.2,1024*0.2), repel =T)+NoLegend()+NoAxes()
DimPlot_idx(reSeuList$seu_strom, group.by = "subtypes00", prefix.index = 72:96,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_strom$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.5,1024*0.5), repel =T)+NoLegend()+NoAxes()
DimPlot_idx(reSeuList$seu_germ, group.by = "subtypes00", prefix.index = 97:105,
            cols = as.vector(mycolortmp[levels(reSeuList$seu_germ$subtypes00)]),
            raster = T, raster.dpi=c(1024*0.2,1024*0.2), repel =T)+NoLegend()+NoAxes()
