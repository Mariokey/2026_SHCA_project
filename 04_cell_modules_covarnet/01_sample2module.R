# ============================================================================
# Sample-to-module assignment and cross-organ module analysis
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Identification of cellular modules"
# Source     : extracted from Analysis240525.r (L1462-1542, L5459-5653)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
sample2module=function(object, vSoi, clusteruse="seurat_clusters", assayuse="SCT", datause="data", result_list=NULL){

DEG_clusterUniq=function(tmpmarkers, minGnum=10){

  final.clusters <- tmpmarkers %>% dplyr::count(cluster) %>% dplyr::filter(n>minGnum) %>% dplyr::pull(cluster)
  markers.final <- tmpmarkers %>% dplyr::filter(cluster %in% final.clusters)
  # Remove duplicated genes
  markers.final$rows <- 1:nrow(markers.final)
  dup <- markers.final[duplicated(markers.final$gene), ]$gene
  remove <- purrr::map(.x=dup, .f=function(i){
    keep <- markers.final %>%
      dplyr::filter(gene== i) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      head(1) %>%
      dplyr::pull(rows)

    remove <- markers.final %>%
      dplyr::filter(gene==i) %>%
      dplyr::filter(rows!=keep) %>%
      dplyr::pull(rows)

    return(remove)
  }) %>% unlist()
  if(length(remove)>=1){
    markers.final <- markers.final[markers.final$rows %!in% remove, ]
  }else{
    markers.final <- markers.final
  }
  final.clusters <- markers.final %>% dplyr::count(cluster) %>% dplyr::filter(n>minGnum) %>% dplyr::pull(cluster)
  markers.final <- markers.final %>% dplyr::filter(cluster %in% final.clusters)
  rownames(markers.final) <- markers.final$gene

  return(markers.final)
}

outRawSeuList=list()
outSeuList=list()
outDegList=list()

for(Soi in vSoi){

  print(paste0(date(), " | Processing ",Soi, " ..."))

  tmpobj=subset(object, sample %in% c(Soi))

  tmpobj=MultiMrun(tmpobj, harmony = 0, mybatch = "sample", sctRun = TRUE)

  # Find markers for banksy clusters (default 0.5 resolution above)
  DefaultAssay(tmpobj) <- assayuse
  Idents(tmpobj) <- clusteruse
  tmpDeg <- FindAllMarkers(tmpobj, test.use = "t", assay = assayuse, slot = datause, logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05, only.pos = TRUE)
  tmpDeg = DEG_clusterUniq(tmpDeg, minGnum = 5)
  tmpDeg$cluster=droplevels(tmpDeg$cluster)

  outRawSeuList[[Soi]]=tmpobj

  select_cluster <- unique(tmpDeg$cluster)
  tmpobj_keep <- subset(tmpobj, seurat_clusters %in% select_cluster)

  outSeuList[[Soi]]=tmpobj_keep
  outDegList[[Soi]]=tmpDeg

  print(paste0(date(), " | ",Soi, " ... done!"))
  if(!is.null(result_list)){
    result_list=list(rawdata=outRawSeuList, data=outSeuList, deg=outDegList)
  }
}

return(list(rawdata=outRawSeuList, data=outSeuList, deg=outDegList))
}

rmv_CTs=c("rmv_BS","rmv_Intestine","rmv_skin","rmv_spleen","rmv_UDT",
          )
subtable=table(inobj$MYL_SubClass04, inobj$Organ)
DEG_clusterUniq=function(tmpmarkers, minGnum=10){

final.clusters <- tmpmarkers %>% dplyr::count(cluster) %>% filter(n>minGnum) %>% pull(cluster)
markers.final <- tmpmarkers %>% filter(cluster %in% final.clusters)
# Remove duplicated genes
markers.final$rows <- 1:nrow(markers.final)
dup <- markers.final[duplicated(markers.final$gene), ]$gene
remove <- purrr::map(.x=dup, .f=function(i){
  keep <- markers.final %>%
    dplyr::filter(gene== i) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    head(1) %>%
    pull(rows)

  remove <- markers.final %>%
    filter(gene==i) %>%
    filter(rows!=keep) %>%
    pull(rows)

  return(remove)
}) %>% unlist()
if(!is.null(remove)){
markers.final <- markers.final[-c(remove), ]
}else{
  markers.final <- markers.final
}
final.clusters <- markers.final %>% dplyr::count(cluster) %>% filter(n>minGnum) %>% pull(cluster)
markers.final <- markers.final %>% filter(cluster %in% final.clusters)
rownames(markers.final) <- markers.final$gene

return(markers.final)
}
outDegList=list()
outSeuList=list()
for(org in colnames(subtable)[colSums(subtable)>0]){
tmp_subdata = subset(inobj, Organ %in% c(org))
tmp_subdata = subset(tmp_subdata, MYL_SubClass04 %!in% rmv_CTs)
tmp_subdata$MYL_SubClass04=droplevels(tmp_subdata$MYL_SubClass04)
tmp_subdata$Organ=droplevels(tmp_subdata$Organ)

Idents(tmp_subdata)="MYL_SubClass04"
tmpDeg <- FindAllMarkers(tmp_subdata,
                       test.use = "t", assay = "RNA", slot = "data", logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05, only.pos = TRUE)
tmpDeg = DEG_clusterUniq(tmpDeg, minGnum = 5)
tmpDeg$cluster=droplevels(tmpDeg$cluster)

select_cluster <- unique(tmpDeg$cluster)
tmpobj_keep <- subset(tmp_subdata, MYL_SubClass04 %in% select_cluster)

outSeuList[[org]]=tmpobj_keep
outDegList[[org]]=tmpDeg
}

seudegList=readRDS("./MP_crossOrgan_module.seudegList_0828.rds")
outSeuList=seudegList$MP_outSeuList
outDegList=seudegList$MP_outDegList

library("tibble")
library("purrr")

keep_n=100
listDE_Banksy_clean=outDegList
all_program_list <-lapply(names(listDE_Banksy_clean),function(i){
  df <- listDE_Banksy_clean[[i]]
  df = df[!grepl("^RP[LS]|^MT-|^ATP",df$gene),]
  df <- df %>% dplyr::group_by(cluster) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::top_n(n = keep_n,wt = avg_log2FC*(1+pct.1-pct.2))  ###p_val_adj
  df <- df %>% dplyr::mutate(Cluster=paste0(i,"_",cluster))
  df_list=tibble_row()
  for(x in unique(df$Cluster)){
    tdf_list=tibble_row()
    tdf_list$Cluster=x
    tdf_list$list=list(df[df$Cluster %in% x,]$gene)
    df_list=rbind(df_list, tdf_list)
  }
  df_list <- df_list %>% dplyr::mutate(len=map_int(df_list$list,function(x){
    length(x)
  }))
  df_list <- df_list %>% dplyr::filter(len>=10)
  return(df_list)
})
all_program_list <- do.call(rbind,all_program_list)
obj_bksy=FBinobj
all_program_list_feed=all_program_list$list
names(all_program_list_feed)=all_program_list$Cluster

obj_bksy=AddModuleScore(obj_bksy, assay = "RNA", name = "PROG.", features = all_program_list_feed, ctrl = 100)
MtxC2Prog=obj_bksy@meta.data[,colnames(obj_bksy@meta.data)[grepl("^PROG\\.", colnames(obj_bksy@meta.data))]]
colnames(MtxC2Prog)=names(all_program_list_feed)
MtxC2Prog=scale(MtxC2Prog)

corMtx=cor(MtxC2Prog)
nModules=8
heattmp=pheatmap(corMtx, cutree_rows = nModules, cutree_cols = nModules,
                 clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 cluster_rows = T, cluster_cols = T, color=viridis::viridis(100), treeheight_row = 0, treeheight_col = 30, name="Cor.")
runFig(paste0(OutputDir,"heatplot.all_FBsubs"), heattmp, 15, 13.5)

rmv_module=c("Module1","Module4") #NULL
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

inobj$Organ=droplevels(inobj$Organ)
keep_CTs=clu_df[clu_df$Cluster %in% unique(clu_df$Cluster),]$ID

for(org in levels(inobj$Organ)){

  xxxtmp=c("Skin_Mph02_STARD13","Skin_LC01_CD207")
  xxxcts= gsub("^[A-Za-z]+_","",xxxtmp) %>% unique()
  tmp_subdata = subset(inobj, Organ %in% c(org))
  tmp_subdata = subset(tmp_subdata, MYL_SubClass04 %in% xxxcts)
  tmp_subdata$Organ=droplevels(tmp_subdata$Organ)
  tmp_subdata$MYL_SubClass04=droplevels(tmp_subdata$MYL_SubClass04)

  Idents(tmp_subdata)="MYL_SubClass04"
  tmpDeg <- FindAllMarkers(tmp_subdata,
                           test.use = "t", assay = "RNA", slot = "data", logfc.threshold = 0.1,
                           min.pct = 0.1, return.thresh = 0.05, only.pos = TRUE)
  tmpDeg = DEG_clusterUniq(tmpDeg, minGnum = 5)
  tmpDeg$cluster=droplevels(tmpDeg$cluster)

  select_cluster <- unique(tmpDeg$cluster)
  tmpobj_keep <- subset(tmp_subdata, MYL_SubClass04 %in% select_cluster)

  outSeuList[[org]]=tmpobj_keep
  outDegList[[org]]=tmpDeg

}

rmv_module=c("Module1","Module4")
MtxC2Prog_sub=scale(MtxC2Prog[,clu_df[clu_df$Cluster %!in% rmv_module,]$ID])
corMtx=cor(MtxC2Prog_sub)
nModules=6
heattmp=pheatmap(corMtx, cutree_rows = nModules, cutree_cols = nModules,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "complete", #"ward.D2",
                 cluster_rows = T, cluster_cols = T, color=viridis::viridis(100), treeheight_row = 10, treeheight_col = 10, name="Cor.")
runFig(paste0(OutputDir,"heatplot.keep_MPsubs"), heattmp, 9.5, 9)
