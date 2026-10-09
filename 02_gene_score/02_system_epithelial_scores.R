# ============================================================================
# Per-system epithelial signature scores and system-level dimension plots
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Gene score analysis"
# Source     : extracted from Analysis240525.r (L10065-10166, L14679-14745)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
#Epi_undef1 as CNTN4
subsysList=readRDS("./subsysList.241220.rds")
tmpseu=subsysList$RS
table(tmpseu$subtypes012)
#"RS"   "US"   "DS"   "SKIN"
golist=list()
golist[["RS"]]=c("GO:1990748", "GO:0009060", "GO:0019731", "GO:0001894", "GO:0060271", "GO:0003341", "GO:0035082", "GO:0042073", "GO:0003351", "GO:0003002", "GO:0060287", "GO:0008544", "GO:0022612", "GO:0060562", "GO:0061138", "GO:0007411", "GO:0097485", "GO:0007409", "GO:0048880", "GO:0060541", "GO:0050808", "GO:0007369")
golist[["US"]]=c("GO:0008544", "GO:0042060", "GO:0022407", "GO:0002065", "GO:0072001", "GO:0001822", "GO:0044782", "GO:0072073", "GO:0097401", "GO:0060993", "GO:0001823", "GO:0072164", "GO:0061326", "GO:0072080", "GO:0001657", "GO:0072017", "GO:0007163", "GO:0050808")
golist[["DS"]]=c("GO:0002181", "GO:0006520", "GO:0010975", "GO:0006892", "GO:0046777", "GO:0007264", "GO:0006631", "GO:0006631", "GO:0006956", "GO:0006066", "GO:0016042", "GO:0006066", "GO:0015986", "GO:1902652", "GO:0008203", "GO:0006869", "GO:0007015", "GO:0015711", "GO:0030216", "GO:0008544", "GO:0050892", "GO:0009913", "GO:0043588", "GO:0006639", "GO:0046486", "GO:0006641", "GO:0022600", "GO:0009225", "GO:0043413", "GO:0050804", "GO:0051047", "GO:0099518", "GO:0046942", "GO:0007586", "GO:0098856", "GO:0015718", "GO:1905039", "GO:0015918", "GO:0030299", "GO:0044241", "GO:0015908", "GO:0042445", "GO:0015850", "GO:0008643", "GO:0034308", "GO:0007586", "GO:0001523", "GO:0019882", "GO:0042060", "GO:0043087", "GO:0022600", "GO:0007163", "GO:0007586", "GO:0001894")
#golist[["DS"]]=names(tmplist) ### 49
golist[["SKIN"]]=c("GO:0043588", "GO:0048732", "GO:0008544", "GO:0022612", "GO:0098773", "GO:0043589", "GO:0030216", "GO:0061436", "GO:0031424", "GO:0042303", "GO:0042633", "GO:0001942")
### Stem;CNTN4;Tuft
#golist[["ADD"]]=names(tmplist) ### 15
### update in 20250219
golist[["ADD"]]=c("GO:0015711", "GO:0010035", "GO:0048762", "GO:0007586", "GO:0030288", "GO:0006887", "GO:0030193", "GO:0030316", "GO:0030363", "GO:0016055", "GO:0001701", "GO:0019827", "GO:0048863", "GO:0030154", "GO:0044763", "GO:0048699", "GO:0030235", "GO:0010550", "GO:0006826", "GO:0072323", "GO:0050905", "GO:0015383", "GO:0015384", "GO:0003927", "GO:0042534", "GO:0046874", "GO:0035437", "GO:0071829", "GO:0030168", "GO:0008233", "GO:0006412", "GO:0010629", "GO:0050905", "GO:0043269", "GO:0033036", "GO:0030054", "GO:0005262", "GO:0007586", "GO:0015801", "GO:0050905", "GO:0015267", "GO:0031053", "GO:0004565", "GO:0060025", "GO:0046863", "GO:0046883", "GO:0030175", "GO:0030855", "GO:0046903", "GO:0031016", "GO:0030295", "GO:0030257", "GO:0031017", "GO:0030429", "GO:0030672", "GO:0001889", "GO:0006584", "GO:0001889", "GO:0031019", "GO:0097367", "GO:0030295", "GO:0035294", "GO:0061493", "GO:0045122", "GO:0035006", "GO:0045785", "GO:0035115", "GO:0061469", "GO:0061467", "GO:0032507", "GO:0032306", "GO:0007600", "GO:0009719", "GO:0061475", "GO:0010629", "GO:0061476", "GO:0051301", "GO:0016477", "GO:0009888", "GO:0030855", "GO:0061478", "GO:0008283", "GO:0030154", "GO:0008544", "GO:0061428", "GO:0061480", "GO:0008283", "GO:0008544", "GO:0045619", "GO:0030198", "GO:0007155", "GO:0061481", "GO:0045619", "GO:0008544", "GO:0008219", "GO:0061428", "GO:0031420", "GO:0061496", "GO:0016337", "GO:0045619", "GO:0008544", "GO:0031420", "GO:0033554", "GO:0040029", "GO:0060496", "GO:0040056", "GO:0040028", "GO:04914", "GO:0060362", "GO:0040029", "GO:0040028", "GO:0032502", "GO:0040056", "GO:04914", "GO:0060362", "GO:0061407", "GO:0061408", "GO:0046903", "GO:0030198", "GO:0030855", "GO:04914", "GO:0042438", "GO:0030316", "GO:0046148", "GO:04914", "GO:0030855", "GO:0016477", "GO:0030855", "GO:0061480", "GO:0008283", "GO:0008544", "GO:0030198", "GO:0030324", "GO:0030154", "GO:0017144", "GO:0097367", "GO:0030324", "GO:0045126", "GO:0030198", "GO:0044782", "GO:0003341", "GO:0030855", "GO:0007163", "GO:0030036", "GO:0042296", "GO:0001822", "GO:0030855", "GO:0072262", "GO:0060362", "GO:0030198", "GO:0035439", "GO:0001822", "GO:0072262", "GO:0030855", "GO:0019829", "GO:0001993", "GO:0006811", "GO:0030855", "GO:0007163", "GO:0030198", "GO:0060362", "GO:0010631", "GO:0048870", "GO:0001541", "GO:0042445", "GO:0030154", "GO:0030198", "GO:0048870", "GO:0048646", "GO:0035976", "GO:0016477", "GO:0007163", "GO:0030198", "GO:0060362", "GO:0002443", "GO:0048666", "GO:0008033", "GO:0030425", "GO:0097458", "GO:0007155", "GO:0007411")

SysID=names(subsysList)
GO_DATA <- clusterProfiler:::get_GO_data(org.Hs.eg.db, "BP", "SYMBOL")

go_vec=unique(golist %>% unlist())
tmplist=list()
for(i in go_vec){
  tmplist[[i]]=GO_DATA$PATHID2EXTID[[i]]
}
for(x in SysID){
  subsysList[[x]]=AddModuleScore(subsysList[[x]], features = tmplist,
                                 ctrl = 100, name = paste0(x,"_GO."))
  print(x)
  print(tail(colnames(subsysList[[x]]@meta.data), n=8))
}

x01=colnames(subsysList$RS@meta.data)[grepl("_GO.", colnames(subsysList$RS@meta.data))]
x01v=names(tmplist)
names(x01v)=x01
x01v_name=GO_DATA$PATHID2NAME[as.vector(x01v)]
x02=colnames(subsysList$US@meta.data)[grepl("_GO.", colnames(subsysList$US@meta.data))]
x02v=names(tmplist)
names(x02v)=x02
x02v_name=GO_DATA$PATHID2NAME[as.vector(x02v)]
x03=colnames(subsysList$DS@meta.data)[grepl("_GO.", colnames(subsysList$DS@meta.data))]
x03v=names(tmplist)
names(x03v)=x03
x03v_name=GO_DATA$PATHID2NAME[as.vector(x03v)]
x04=colnames(subsysList$SKIN@meta.data)[grepl("_GO.", colnames(subsysList$SKIN@meta.data))]
x04v=names(tmplist)
names(x04v)=x04
x04v_name=GO_DATA$PATHID2NAME[as.vector(x04v)]

subsysList=readRDS("./subsysListNEW.250219.rds")
sys_v=names(subsysList)
Sys_plot_list=list()
legend_df_vec=c()
pdf("./sys_dimplot.pdf", width = 6.5, height = 4.5)
for(sv in sys_v){
    tmpd=subsysList[[sv]]

prefix_id=names(CT107_color)[32:61]
prefix_id_order=1:length(prefix_id)
names(prefix_id_order)=prefix_id
prefix_id=prefix_id_order[levels(tmpd$refine_subtypes01)]
legend_df_vec=c(legend_df_vec, prefix_id)
p=scPioneer::DimPlot_idx(tmpd, group.by = "refine_subtypes01", cols=as.vector(CT107_color[levels(tmpd$refine_subtypes01)]),
                         raster=TRUE, raster.dpi=c(350,350),
                         prefix.index=as.vector(prefix_id))+NoAxes() #+NoLegend()
print(p+ggtitle(label = sv))

}
newdf=data.frame(celltype=names(legend_df_vec), index=as.vector(legend_df_vec)) %>% unique()
newdf=newdf[order(newdf$index),]
newdf$x=c(rep(1,8),rep(2,12),rep(3,8))
newdf$y=c(10:3, 12:1, 10:3)
newdf$annotation.idx=paste0(newdf$index," - ", newdf$celltype)
ggplot(newdf)+
  geom_point(aes(x=x-0.1, y=y, fill=celltype, color=celltype), #color="white",
             shape=21,
             size=5)+
  scale_fill_manual(values = CT107_color)+
  scale_color_manual(values = CT107_color)+
  geom_text(aes(x=x, y=y, label=annotation.idx ), color="black", hjust=0)+  ###color=annotation
  theme_nothing()+xlim(c(0.5,5))
dev.off()
