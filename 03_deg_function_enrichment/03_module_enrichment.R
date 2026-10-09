# ============================================================================
# Functional enrichment of the cross-organ cellular modules
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Function annotation"
# Source     : extracted from Analysis240525.r (L5655-5835)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
rmv_module=NULL
rcl.list=row_order(draw(heattmp)) %>% as.list()
clu_df = data.frame()
clu_df = lapply(1:length(rcl.list), function(i){
  idid=rownames(corMtx)[rcl.list[[i]]] %>% as.vector()
  out <- data.frame(ID = idid, Cluster = rep(paste0("Module", i), length(idid)))
  return(out)
}) %>% do.call(rbind, .) %>% as.data.frame()
clu_df$IDsplit=clu_df$ID
clu_df$IDsplit=gsub("_\\S+", "", clu_df$IDsplit)
stat_clu=table(clu_df$Cluster, clu_df$IDsplit) %>% as.matrix()
clu_df=clu_df[clu_df$Cluster %in% rownames(stat_clu[rowSums(stat_clu>0)>0,]),]
clu_df=clu_df[clu_df$Cluster %!in% rmv_module,]

enrich_gene_list=list()
ModuleFun_raw = lapply(unique(clu_df$Cluster), function(i, alist=enrich_gene_list){
  df_list=tibble_row()
  IDvec=clu_df[clu_df$Cluster %in% i,]$ID
  IDvec_gene=all_program_list_feed[IDvec] %>% unlist() %>% as.vector() %>% table() %>% sort()
  tmpDFdata=data.frame(x=1:length(IDvec_gene), gene=names(IDvec_gene) %>% rev(), counts=rev(IDvec_gene %>% as.vector()))
  alist[[i]]=ggplot(tmpDFdata)+geom_point(aes(x,counts), color="lightgrey")+
    ggrepel::geom_text_repel(data=tmpDFdata[tmpDFdata$counts>quantile(tmpDFdata$counts, probs=0.9),],
                             aes(x,counts, label=gene), force_pull = 1, force = 5,direction = "y",
                             nudge_x = 250, min.segment.length = 0, max.overlaps = 8, nudge_y = -5)+theme_cowplot()+xlab("")
  if(max(IDvec_gene)>2){
    IDvec_gene=names(IDvec_gene[IDvec_gene>=2])
  }else{
    IDvec_gene=names(IDvec_gene)
  }

  tdf_list=tibble_row()
  tdf_list$ModuleID=i
  tdf_list$list=list(IDvec_gene)
  tdf_list$funlist=list(RunHsGO(IDvec_gene, trans = TRUE))

  df_list=rbind(df_list, tdf_list)

  return(list(df_list, alist))
})
ModuleFun=ModuleFun_raw[[1]]
enrich_gene_list=ModuleFun_raw[[2]]
ModuleFun=rbind(
  ModuleFun_raw[[1]][[1]],
  ModuleFun_raw[[2]][[1]],
  ModuleFun_raw[[3]][[1]],
  ModuleFun_raw[[4]][[1]],
  ModuleFun_raw[[5]][[1]],
  ModuleFun_raw[[6]][[1]],
  ModuleFun_raw[[7]][[1]],
  ModuleFun_raw[[8]][[1]]
          )

godf=data.frame()
for(x in 1:dim(ModuleFun)[1]){
  mmmdf=as.data.frame(ModuleFun[x,]$funlist[[1]]$d)
  mmmdf$Module=rep(ModuleFun[x,]$ModuleID, dim(mmmdf)[1])
  godf=rbind(godf,mmmdf)
}
godf=godf %>% group_by(Module) %>% top_n(n = 8, wt = -pvalue)
godf$Module=factor(godf$Module, levels = paste("Module",1:length(unique(godf$Module)), sep = ""))
godf$Description=factor(godf$Description, levels = unique(godf$Description))
ggplot(godf, aes(Module, Description))+
  geom_point(aes(size=Factor,fill=-log10(pvalue)), shape=22)+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = 45)+
  xlab("")+ylab("")+
  theme_cowplot()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.y = element_blank()) #+

write.xlsx(godf, paste0(OutputDir, "/MPmajor_5Module.FuncEnrich.xlsx"))
write.xlsx(godf, paste0(OutputDir, "/MPresident_6Module.FuncEnrich.xlsx"))

write.xlsx(obj_bksy_degs, paste0(OutputDir, "/MPresident_6Module.byDEGs.xlsx"))
write.xlsx(godf, paste0(OutputDir, "/MPresident_6Module.byDEGs2Func.xlsx"))

write.xlsx(godf, paste0(OutputDir, "/FB_8Module.FuncEnrich.xlsx"))
fboverlapGenes_list=list()
for(x in 1:8){
  fboverlapGenes_list[[x]]=ModuleFun_raw[[x]][[2]][[1]]
}
xxx=plot_grid(plotlist = fboverlapGenes_list, ncol = 2, labels = paste0("Module",1:8))
runFig(paste0(OutputDir,"intersectGenes_FBsubs_8Modules"), xxx, 10, 18)

OutputDir="./AnaRes20240829/"
dir.create(OutputDir, recursive = T)

obj_bksy$organ_ct04=paste0(obj_bksy$Organ,"_",obj_bksy$MYL_SubClass04)
Module_list=list()
for(MMM in unique(clu_df$Cluster)){
Module_list[[MMM]]=clu_df[clu_df$Cluster %in% MMM,]$ID   ###
}
negVec=NULL
obj_bksy$organ_ct04=factor(obj_bksy$organ_ct04, levels = unique(obj_bksy$organ_ct04))
for(x in 1:length(Module_list)){
  obj_bksy=modifyAno(obj_bksy, anoCol = "organ_ct04",
                  cellid = setdiff(subset(obj_bksy, organ_ct04 %in% Module_list[[x]])$cellid, negVec),
                  newID = names(Module_list[x]))
}
unique(obj_bksy$organ_ct04)

obj_bksy$organ_ct04_bak=paste0(obj_bksy$Organ,"_",obj_bksy$MYL_SubClass04)
obj_bksy$organ_ct04_bak=factor(obj_bksy$organ_ct04_bak, levels = unique(obj_bksy$organ_ct04_bak))
obj_bksy_sub=subset(obj_bksy, organ_ct04 %in% c("Module1","Module2","Module3","Module4","Module5","Module6"))
obj_bksy_sub=subset(obj_bksy_sub, organ_ct04_bak %!in% c("Testis_MphX_OvaryTestis.C1QB",
                                                        "Ovary_MphX_OvaryTestis.C1QB",
                                                        "Airway_MphX_Airway.LPL",
                                                        "Lung_Mph01_AGRP",
                                                        "FU_DC01_CD1C",
                                                        "UB_MphX_UB.C1QB"))

obj_bksy_sub$organ_ct04=droplevels(obj_bksy_sub$organ_ct04)
Idents(obj_bksy_sub)="organ_ct04"
obj_bksy_degs=FindAllMarkers(obj_bksy_sub,
                             group.by = "organ_ct04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
DotPlot(obj_bksy_sub, features = (obj_bksy_degs %>% group_by(cluster) %>% top_n(8, wt =avg_log2FC))$gene %>% unique(),
        group.by = "organ_ct04") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")

obj_bksy_sub <- ScaleData(obj_bksy_sub, features = rownames(obj_bksy_sub))
obj_bksy_sub$cellnames=rownames(obj_bksy_sub@meta.data)
cellID <- subsetID2(obj_bksy_sub@meta.data[,c('organ_ct04','cellnames')],expected.cell = 400, seed=1024)
subobj <- subset(obj_bksy_sub, cells = cellID$Cellnames)
p1 <- DoHeatmap(subobj,  #obj_bksy, #slot="data",
                features = (obj_bksy_degs %>% group_by(cluster) %>% top_n(14, wt =avg_log2FC*(0.5+pct.1-pct.2)))$gene %>% unique(),
                size = 3,angle = 30,group.bar.height = 0.05)
p1[[1]]$layers[[2]]$show.legend=F
p=p1[[1]]+scale_fill_gradient2(low = "#141464", mid = "white", high = "#78050F", midpoint = 0, limits=c(-2.5,2.5))

obj_bksy_funcs=enrichFUNgo(ud_table = obj_bksy_degs,intermode = FALSE,enriches = "GO")
godf=obj_bksy_funcs$dflist$GO %>% group_by(celltype) %>% top_n(n = 100, wt = -pvalue)
godf=obj_bksy_funcs$dflist$GO %>% group_by(celltype) %>% top_n(n = 10, wt = -pvalue)
godf$Module=godf$celltype
godf$Module=factor(godf$Module, levels = paste("Module",1:length(unique(godf$Module)), sep = ""))
godf$Description=factor(godf$Description, levels = unique(godf$Description))
ggplot(godf, aes(Module, Description))+
  geom_point(aes(size=Factor,fill=-log10(pvalue)), shape=22)+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = 45)+
  xlab("")+ylab("")+
  theme_cowplot()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.y = element_blank()) #+

saveRDS(obj_bksy, file="./AnaRes20240913/MPsub_only.obj_bksy.rds")
###
###
### %>% gsub("\S+_","",.)
###
OutputDir="./AnaRes20240913/"
dir.create(OutputDir, recursive = T)

rename_list=list()
rename_list[["Mph01_STARD13"]]=clu_df[clu_df$Cluster %in% "Module1",]$ID
rename_list[["Mph02_SELENOP"]]=clu_df[clu_df$Cluster %in% "Module2",]$ID
rename_list[["Mph03_NR4A1"]]=clu_df[clu_df$Cluster %in% "Module3",]$ID   ### MoMP
rename_list[["Mph04_MCEMP1"]]=clu_df[clu_df$Cluster %in% "Module4",]$ID
rename_list[["MphX_VDFU.C1QB"]]=clu_df[clu_df$Cluster %in% "Module5",]$ID
rename_list[["Mono01_CD14"]]=c()
rename_list[["DC01_CD1C"]]=c()

inobj$MYL_SubClass05=inobj$MYL_SubClass04
negVec=NULL
for(x in 1:length(rename_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass05",
                  cellid = setdiff(subset(inobj, MYL_SubClass05 %in% rename_list[[x]])$cellid, negVec),
                  newID = names(rename_list[x]))
}
