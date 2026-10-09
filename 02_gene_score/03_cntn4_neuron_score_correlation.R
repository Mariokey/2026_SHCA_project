# ============================================================================
# Correlation between CNTN4+ cell signature score and neuron-receptor score
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Gene score analysis"
# Source     : extracted from Analysis240525.r (L14900-15008)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
tmpDF=colPercents(tmpxxx)[1:dim(tmpxxx)[1],]
color_state=CT107_color
p=ggplot(tmpDF %>% melt() %>% set_colnames(.,c("Celltype","Group","Fraction")))+
  geom_bar(aes(Group, Fraction, fill=as.factor(Celltype), color=as.factor(Celltype)), stat = "identity", position = position_stack())+
  scale_fill_manual(values = color_state)+scale_color_manual(values = color_state)+
  xlab("")+theme_minimal()+myfont+NoLegend()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
runFig(paste0(OutputDir,"/", "Fig2.episub_dodgeBar_withRPS"), p, 12.5, 2.15)

xxxdata=subset(Mobj_keep, MajorClass02 %in% c("Epithelial Cells","Endothelial Cells"))
xxxdata$subtypes015=droplevels(xxxdata$subtypes015)
pdf("./NE_syn_degrade_AvgScore_splitby_LocCt.pdf", width = 18, height = 15)
  tmp_df1=xxxdata@meta.data[,c("nesyn1","Location_short","subtypes015")] %>%
    dcast(Location_short~subtypes015, value.var = "nesyn1", fun=mean) #tidyr::pivot_wider()
  tmp_df11=tmp_df1[,-1]
  rownames(tmp_df11)=tmp_df1[,1]
  tmp_df11=as.matrix(tmp_df11)
  tmp_df11[is.na(tmp_df11)]=min(tmp_df11,na.rm = TRUE)
  p=pheatmap(tmp_df11, name=paste0("avg. score\n","NE_syn&deg"),
             cluster_rows = F, cluster_cols = F)
  draw(p)
dev.off()

afterQC_df=read.xlsx("./data/upload_oss/uploadRes_240425/All_metadata_91w.0426.xlsx")
afterQC_df00=afterQC_df[,c("Sample_ID","sample_old","Location_short","n_genes_by_counts","total_counts","total_counts_mt","pct_counts_mt","log10genes_by_umi","n_counts","n_genes")] %>% mutate(dataset=rep("after_QC",dim(afterQC_df)[1]))
beforeQC_df=srt_meta
beforeQC_df00=left_join(beforeQC_df, unique(afterQC_df[,c("Sample_ID","sample_old","Location_short")]), by=c("sample"="sample_old"))
colnames(beforeQC_df00)=gsub("^sample$","sample_old",colnames(beforeQC_df00))
beforeQC_df00=beforeQC_df00[,c("Sample_ID","sample_old","Location_short","n_genes_by_counts","total_counts","total_counts_mt","pct_counts_mt","log10genes_by_umi","n_counts","n_genes")] %>% mutate(dataset=rep("before_QC",dim(beforeQC_df00)[1]))
joinQC_df=rbind(beforeQC_df00,afterQC_df00)
sample_level=unique(Mobj_keep@meta.data[,c("Location_short","Sample_ID")]) %>% arrange_at("Location_short")
sample_level=sample_level$Sample_ID %>% as.vector()
location_level=levels(Mobj_keep$Location_short)
joinQC_df$Sample_ID=factor(joinQC_df$Sample_ID, levels = rev(sample_level))
joinQC_df$Location_short=factor(joinQC_df$Location_short, levels = location_level)
joinQC_df$dataset=factor(joinQC_df$dataset, levels = c("before_QC","after_QC"))
mycolumn="n_genes"
mycolumn="n_counts"
mycolumn="pct_counts_mt"
tmpP=ggplot(joinQC_df)+
  geom_boxplot(aes_string("Sample_ID", mycolumn, fill="dataset"), outlier.colour = "white", outlier.size = 0, outlier.alpha = 0)+
  scale_fill_manual(values=c("before_QC"="grey","after_QC"="red"))+scale_color_manual(values=c("before_QC"="grey","after_QC"="red"))+
  theme_minimal_grid()+coord_flip()+xlab("")#+ylim(c(0,50000))
runFig(paste0(FinalFig,"/", mycolumn, "_QCboxplot"), tmpP, 8, 24)

Marker_Lev1=c("KRT1","KRTDAP","DMKN","EPCAM",
              "PTPRC","CD3E","CD68","CD79A",
              "ACTA2","MYH11","DCN","PDGFRB",
              "PLVAP","CDH5","PECAM1","VWF",
              "TNP1","PRM1","PRM3","DDX4",
              "MKI67","TOP2A","UBE2C","EZH2")
scPioneer::FeaturePlot2(Mobj_keep, features = Marker_Lev1, reduction = "umap",label = F, pt.size = 1, ncol = 4)
ggsave(filename = paste0(FinalFig,"/Mobj_Marker.pdf"),width = 12, height = 16)

  x="Stromal Cells"
  lollipop_plot(
    tmpprop_filter[tmpprop_filter$Cell_class %in% x,] %>% unique(),   #### ==> line 9414
    "Location","logFC",group = "System",
  )#+theme_minimal_vgrid()+NoLegend()
  ggsave(filename = paste0(FinalFig,"/Mobj_CellClass_LocEnrich_",gsub("\\.","",gsub("\\s+\\S+","",x)),".pdf"),width = 4.5, height = 6.5)
#}
    lollipop_plot(
      tmpprop_filter00[tmpprop_filter00$Cell_class %in% x,] %>% unique(),
      "System","logFC",
      dot_color = newColor
    )+theme(axis.line.y = element_line(colour = "grey"))+geom_vline(xintercept = 0, color = "lightgrey",linewidth=0.6)+geom_vline(xintercept = 1) #+NoAxes()
    #+theme_minimal_vgrid()+NoLegend()
    ggsave(filename = paste0(FinalFig,"/Mobj_SysEnrich_",gsub("\\.","",gsub("\\s+\\S+","",x)),".OnlySys.pdf"),width = 2.75, height = 2.75)
 # }
#}

tmpP=ggplot(CNTN4_data_ori@meta.data[,c("neuron1","CNTN4_2score1","types")])+geom_point(aes(neuron1, CNTN4_2score1, color=types), alpha=0.3)
tmpP_splitplot=list()
for(x in names(table(tmpP$data$types))){
  tmpP_data=tmpP$data[tmpP$data$types %in% c(x),]
  print(x)
  y=cor.test(tmpP_data$neuron1, tmpP_data$CNTN4_2score1)
  yp=round(y$p.value, digits = 3)
  yr=round(y$estimate, digits = 3)
  tmpP_splitplot[[x]]=ggplot(tmpP_data)+
    geom_point(aes(neuron1, CNTN4_2score1), shape=21, size=1, color="darkgrey", alpha=1)+theme_cowplot()+
    geom_smooth(aes(neuron1, CNTN4_2score1), method="lm")+
    ggtitle(paste0(x,"; cor: ",yr,"; p-value:",yp))+
    xlab("Neuron receptor score")+ylab("CNTN4_Cell signature score")
}
plot_grid(plotlist = tmpP_splitplot,  ncol = 1, axis="tblr", align = "hv")
ggsave(filename = paste0(FinalFig,"/CNTN4_NeuronR_Cor.pdf"),width = 5.5, height = 25)
tmpP+theme_cowplot()+geom_smooth(aes(neuron1, CNTN4_2score1), method="lm")+
  ggtitle(paste0("cor: ","0.264","; p-value:","2.2e-16"))+
  xlab("Neuron receptor score")+ylab("CNTN4_Cell signature score")
ggsave(filename = paste0(FinalFig,"/CNTN4_NeuronR_AllCombinCor.pdf"),width = 5, height = 4)
