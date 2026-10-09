# ============================================================================
# Epithelial label curation and system-level regrouping
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L3832-4201)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### OBSOLATE/obsolate

###
seu_epi_all=readRDS("./test.epi_seuratAnno.v0525.rds")
###
seu_epi_all$subtypes011=seu_epi_all$SubClass02
seu_epi_all=subset(seu_epi_all, subtypes011 %!in% c("Glia_Cells","Immune_Cells","Neurons","RPS_CS_Epithelial_Cells_2"))
seu_epi_all$subtypes011=droplevels(seu_epi_all$subtypes011)
seu_epi_all=modifyAno(seu_epi_all, anoCol = "subtypes011", oldID = "DS_Epithelial_Cells_1", newID = "Cholangiocytes",
               newOrder = gsub("DS_Epithelial_Cells_1","Cholangiocytes",levels(seu_epi_all$subtypes011)))
seu_epi_all=modifyAno(seu_epi_all, anoCol = "subtypes011", oldID = "RPS_Epithelial_Cells_2", newID = "Granulosa",
               newOrder = gsub("RPS_Epithelial_Cells_2","Granulosa",levels(seu_epi_all$subtypes011)))
seu_epi_all=modifyAno(seu_epi_all, anoCol = "subtypes011", oldID = "RPS_CS_Epithelial_Cells_1", newID = "Epi_undef1",
               newOrder = gsub("RPS_CS_Epithelial_Cells_1","Epi_undef1",levels(seu_epi_all$subtypes011)))
seu_epi_all=modifyAno(seu_epi_all, anoCol = "subtypes011", oldID = "MIX_Epithelial_Cells_1", newID = "Epi_undef2",
               newOrder = gsub("MIX_Epithelial_Cells_1","Epi_undef2",levels(seu_epi_all$subtypes011)))
seu_epi_all$subtypes012=seu_epi_all$subtypes011
seu_epi_all$subtypes012=factor(seu_epi_all$subtypes012, levels =
                                 intersect(Mobj$subtypes011 %>% levels(), seu_epi_all$subtypes011))

resDir="./AnaRes20240704/"
dir.create(resDir, recursive = T)
SYSoi=c("RS", "US", "DS", "SKIN")
minRate=c(0.02, 0.02, 0.02, 0.02) #SKIN:0.01?
plotlist=list()
DEG2FUN_subcts_list=list()
df_DEGs=data.frame()
df_GOs=data.frame()
df_CellNums=data.frame()
for(i in 1:length(SYSoi)){
  #i=3
  mysys=SYSoi[i]

typecol="subtypes012" #paste0(shortCToi[i],"subtype00")
tmpseu=subsysList[[mysys]]
DefaultAssay(tmpseu)="RNA"
Idents(tmpseu)=typecol

deg000=DEG2FUN_subcts_list[[mysys]]$RawDEG[
  abs(DEG2FUN_subcts_list[[mysys]]$RawDEG$avg_log2FC)>0.25 &
  DEG2FUN_subcts_list[[mysys]]$RawDEG$p_val_adj<0.01,]
fun000=DEG2FUN_subcts_list[[mysys]]$FUNdf$GO

deg000$SYS=mysys
df_DEGs=rbind(df_DEGs, deg000)
fun000$SYS=mysys
df_GOs=rbind(df_GOs, fun000)

tmpseu$Location_short=droplevels(tmpseu$Location_short)
cnum000=table(tmpseu$subtypes012, tmpseu$Location_short) %>% as.data.frame()
cnum000$SYS=mysys
df_CellNums=rbind(df_CellNums, cnum000)

}

colorlist=list()
color_range=c(0, 0.125, 0.35, 0.8)
plotlist=list()
for(i in 1:length(SYSoi[c(3,4,1,2)])){
  mysys=SYSoi[c(3,4,1,2)][i]
  myCT=subsysList[[mysys]]$subtypes012 %>% levels()
  tmpColor=rainbow(start = color_range[i],end = color_range[i]+0.15, n = length(myCT))
  names(tmpColor)=myCT
  colorlist[[mysys]]=tmpColor
  plotlist[[mysys]]=DimPlot(subsysList[[mysys]], cols = tmpColor, group.by = "subtypes012", reduction = "umap", label = T, repel = T)+NoLegend()+ggtitle("")
}
plot_grid(plotlist = plotlist, nrow = 1)

df_DEGs_nums=table(df_DEGs$SYS,df_DEGs$celltype) %>% as.data.frame()
df_DEGs_nums$SysCT=paste(df_DEGs_nums$Var1,df_DEGs_nums$Var2, sep = ".")
df_DEGs_nums=df_DEGs_nums[df_DEGs_nums$Freq>0,]
df_DEGs_nums$SysCT=factor(df_DEGs_nums$SysCT, levels=names(unlist(colorlist)))
ggplot(df_DEGs_nums)+geom_bar(aes(SysCT, Freq/1000, fill=SysCT),color="white", stat = "identity")+scale_fill_manual(values = unlist(colorlist))+theme_minimal()+NoLegend()+xlab("")+ylab("#DEGs (x1000)")+theme(axis.text.x = element_blank())

df_CellNums$SYS=factor(df_CellNums$SYS, levels = SYSoi[c(3,4,1,2)])
df_CellNums$Var1=factor(df_CellNums$Var1, levels = levels(seu_epi_all$subtypes012))
df_CellNums$Var2=factor(df_CellNums$Var2, levels = levels(seu_epi_all$Location_short))
ggplot(df_CellNums[df_CellNums$Freq>0,])+geom_bar(aes(Var1, Freq/10000, fill=Var2), color="white", stat = "identity", position = position_stack())+facet_wrap(~SYS, scales = "free_x", nrow = 1)+scale_x_discrete(drop=TRUE)+theme(axis.text.x = element_text(angle = 90))

###
indir="./AnaRes20240704/"
xlsx_files=c("DS_subcts.results.081238.GO.xlsx",
             "SKIN_subcts.results.084829.GO.xlsx",
             "RS_subcts.results.072934.GO.xlsx",
             "US_subcts.results.073759.GO.xlsx")
syss=c("DS","SKIN","RS","US")
df_GOs=data.frame()
for(x in 1:4){
  fun000=read.xlsx(paste0(indir, "/",xlsx_files[x]))
  head(fun000)
  fun000$SYS=syss[x]
  df_GOs=rbind(df_GOs, fun000)
}
df_GOs$GOCT=paste0(df_GOs$ID,"<>",df_GOs$celltype)
df_GOs$SysCT=paste(df_GOs$SYS,df_GOs$celltype, sep = ".")

selected_GOs=list()
selected_GOs[["DS"]]=text2vector('GO:0006520<>Hepatocytes
GO:0006631<>Hepatocytes
GO:0016042<>Hepatocytes
GO:0006869<>Enterocytes
GO:0015711<>Enterocytes
GO:0050892<>Enterocytes
GO:0022600<>Enterocytes
GO:0046942<>Enterocytes
GO:0007586<>Enterocytes
GO:0098856<>Enterocytes
GO:0015718<>Enterocytes
GO:1905039<>Enterocytes
GO:0015918<>Enterocytes
GO:0030299<>Enterocytes
GO:0044241<>Enterocytes
GO:0015908<>Enterocytes
GO:0015850<>Enterocytes
GO:0008643<>Enterocytes
GO:0007586<>Epi_undef2
GO:0022600<>BEST4_Epithelial_Cells
GO:0007163<>Prezymogenic_Cells
GO:0007586<>BEST4_Epithelial_Cells')
selected_GOs[["SKIN"]]=text2vector('GO:0043588<>Basal_Cells
GO:0048732<>Basal_Cells
GO:0008544<>Basal_Cells
GO:0022612<>Basal_Cells
GO:0098773<>Basal_Cells
GO:0043589<>Basal_Cells
GO:0030216<>Early_Suprabasal_Cells
GO:0061436<>Early_Suprabasal_Cells
GO:0031424<>Early_Suprabasal_Cells
GO:0042303<>Vellus_Hair_Follicle
GO:0042633<>Vellus_Hair_Follicle
GO:0001942<>Vellus_Hair_Follicle')
selected_GOs[["RS"]]=text2vector('GO:0060271<>Multiciliated_Cells
GO:0003341<>Multiciliated_Cells
GO:0035082<>Multiciliated_Cells
GO:0042073<>Multiciliated_Cells
GO:0003351<>Multiciliated_Cells
GO:0003002<>Multiciliated_Cells
GO:0060287<>Multiciliated_Cells
GO:0008544<>Papillary_Epithelial_Cells
GO:0022612<>Papillary_Epithelial_Cells
GO:0060562<>Papillary_Epithelial_Cells
GO:0061138<>Papillary_Epithelial_Cells
GO:0007411<>Papillary_Epithelial_Cells
GO:0097485<>Papillary_Epithelial_Cells
GO:0007409<>Papillary_Epithelial_Cells
GO:0048880<>Papillary_Epithelial_Cells
GO:0060541<>Papillary_Epithelial_Cells
GO:0050808<>Epi_undef1
GO:0007369<>Epi_undef1')
selected_GOs[["US"]]=text2vector('GO:0002065<>Papillary_Epithelial_Cells
GO:0072001<>Cortical_Collecting_Duct_Principal_Cells
GO:0072073<>Cortical_Collecting_Duct_Principal_Cells
GO:0097401<>Cortical_Collecting_Duct_Principal_Cells
GO:0060993<>Cortical_Collecting_Duct_Principal_Cells
GO:0001823<>Cortical_Collecting_Duct_Principal_Cells
GO:0072164<>Cortical_Collecting_Duct_Principal_Cells
GO:0061326<>Cortical_Collecting_Duct_Principal_Cells
GO:0072080<>Cortical_Collecting_Duct_Principal_Cells
GO:0001657<>Cortical_Collecting_Duct_Principal_Cells
GO:0072017<>Cortical_Collecting_Duct_Principal_Cells
GO:0007163<>Epi_undef1
GO:0050808<>Epi_undef1')
df_GOs=df_GOs[df_GOs$Description %!in%  c("cytoplasmic translation"),]
df_GOs=df_GOs[!grepl("ribosome|translation",df_GOs$Description),]
df_GOs_top1=df_GOs %>% group_by(SysCT) %>% top_n(n = 1, wt=Factor*(-1*log10(p.adjust)))
df_GOs_top3=df_GOs[df_GOs$GOCT %in% (unlist(selected_GOs) %>% as.vector()),]
myterms=table(df_GOs_top3$Description)#[table(df_GOs_top3$Description)>=2 & table(df_GOs_top3$Description)<=5]
df_GOs_plotdf=df_GOs_top3[df_GOs_top3$Description %in% names(myterms),]
df_GOs_plotdf=rbind(df_GOs_plotdf, df_GOs_top1) %>% unique()
df_GOs_plotdf_matrix=df_GOs_plotdf %>% tidyr::pivot_wider(names_from = "Description", id_cols = "SysCT", values_from = "Factor")
rownnn=df_GOs_plotdf_matrix$SysCT
df_GOs_plotdf_matrix=df_GOs_plotdf_matrix[,2:79] %>% as.matrix()
rownames(df_GOs_plotdf_matrix)=rownnn
df_GOs_plotdf_matrix[is.na(df_GOs_plotdf_matrix)]=0
df_GOs_plotdf_matrix=df_GOs_plotdf_matrix[intersect(names(unlist(colorlist)),rownames(df_GOs_plotdf_matrix)),]
dispMatrix=t(df_GOs_plotdf_matrix)
dispMatrix=round(dispMatrix,digits = 3)
dispMatrix[dispMatrix==0]=""
pheatmap(t(df_GOs_plotdf_matrix), name = "Factor", show_colnames = F, show_rownames = T,
         cluster_rows = F, cluster_cols = F,number_color ="black",
         color = viridis_pal()(100), border_color = NA, number_format = "%.2f")

matches_freq=-log10(tmpdf$p.adjust)
names(matches_freq)=tmpdf$Description
wc_color=rep("black",length(matches_freq))
wordcloud(names(matches_freq), matches_freq, scale=c(3.5, 0.25), min.freq=1, random.order=FALSE, fixed.asp=F, rot.per=0,
          colors=wc_color, ordered.colors=T,
          random.color=F)

library(wordcloud)
