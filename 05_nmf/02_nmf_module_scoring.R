# ============================================================================
# NMF program scoring, correlation and Jaccard clustering of programs
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "NMF analysis"
# Source     : extracted from Analysis240525.r (L13590-13983)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
cm_network <- function (NMFres, Corres, ...) {
  var_args <- list(...)
  top_n <- if (!is.null(var_args[["top_n"]]))
    var_args[["top_n"]]
  else 10
  corr <- if (!is.null(var_args[["corr"]]))
    var_args[["corr"]]
  else 0.2
  fdr <- if (!is.null(var_args[["pval_fdr"]]))
    var_args[["pval_fdr"]]
  else 0.05
  listA = Corres[, c("subCluster1", "majorCluster1")]
  listB = Corres[, c("subCluster2", "majorCluster2")]
  colnames(listA) = c("subCluster", "majorCluster")
  colnames(listB) = c("subCluster", "majorCluster")
  ann = rbind(listA, listB)
  ann = ann[!duplicated(ann$subCluster), ]
  rownames(ann) = ann$subCluster
  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
  }
  sorted_colnames <- order(colnames(w))
  w <- w[, sorted_colnames]
  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster", "cm", "weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>%
    ungroup()
  weight <- weight %>% group_by(cm) %>% arrange(desc(weight)) %>%
    clusterProfiler::slice(1:top_n) %>% ungroup()
  meta_sc = data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                       majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]
  rownames(meta_sc) <- meta_sc$subCluster
  meta_sc <- meta_sc[order(meta_sc$subCluster), ]
  cdt1 <- Corres$correlation > corr
  cdt2 <- Corres$pval_fdr <= fdr
  pl_df <- Corres[cdt1 & cdt2, ]
  n_all <- nrow(meta_sc)
  spe_cutoff <- 1 - ((top_n - 1) * 2 - 1)/((n_all - 1) * 2 -
                                             1)
  cdt3 <- pl_df$spe >= spe_cutoff
  pl_df <- pl_df[cdt3, ]
  graph <- graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- V(graph)$name
  node_global <- data.frame(subCluster = node_global, majorCluster = meta_sc$majorCluster[match(node_global,
                                                                                                meta_sc$subCluster)])
  node_global <- node_global[order(node_global$subCluster),
  ]
  node_each <- matrix(NA, nrow = 0, ncol = 2) %>% as.data.frame()
  colnames(node_each) <- c("cm", "subCluster")
  edge_each <- matrix(NA, nrow = 0, ncol = ncol(pl_df) + 1) %>%
    as.data.frame()
  colnames(edge_each) <- c("cm", colnames(pl_df))
  for (cm in unique(weight$cm)) {
    sub_graph <- subgraph(graph, vids = intersect(V(graph)$name,
                                                  weight$subCluster[weight$cm == cm]))
    sub_graph <- delete.vertices(sub_graph, V(sub_graph)[degree(sub_graph) ==
                                                           0])
    if (length(V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                      (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm,
                                               tmp2))
    }
  }
  node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster,
                                                       meta_sc$subCluster)]
  rownames(edge_each) <- NULL
  weight = merge(node_each, weight)
  weight = weight[order(weight$cm, -weight$weight), ]
  rownames(weight) = as.list(1:nrow(weight))
  cm_network <- list(list(node_global, pl_df), list(node_each,
                                                    edge_each), NMFres, weight, ann)
  names(cm_network) <- c("global", "each", "raw", "filter",
                         "ann")
  names(cm_network$global) <- c("node", "edge")
  names(cm_network$each) <- c("node", "edge")
  return(cm_network)
}

Xgr.igraph_each=function (each, ...) {
  var_args <- list(...)
  Layout <- if (!is.null(var_args[["Layout"]]))
    var_args[["Layout"]]
  else layout_in_circle
  node = each$node
  edge = each$edge
  edge = edge[, c("subCluster1", "subCluster2", "cm", "correlation",
                  "pval", "pval_fdr", "spe", "majorCluster1", "majorCluster2")]
  graph <- graph_from_data_frame(edge, directed = FALSE)
  meta_sc = data.frame(subCluster = c(edge$subCluster1, edge$subCluster2),
                       majorCluster = c(edge$majorCluster1, edge$majorCluster2))
  meta_sc = meta_sc[!duplicated(meta_sc$subCluster), ]
  rownames(meta_sc) = meta_sc$subCluster
  V(graph)$majorCluster <- meta_sc[V(graph)$name, "majorCluster"]
  colors <- RColorBrewer::brewer.pal(length(unique(V(graph)$majorCluster)),
                                     "Set3")
  names(colors) = unique(V(graph)$majorCluster)
  V(graph)$frame.color <- V(graph)$color <- colors[as.vector(V(graph)$majorCluster)]
  col_fun <- circlize::colorRamp2(breaks = quantile(E(graph)$spe,
                                                    probs = seq(0, 1, length = 6)), colors = (ggsci::pal_material(palette = "grey",
                                                                                                                  n = 10))(10)[3:8])
  E(graph)$color <- col_fun(E(graph)$spe)
  E(graph)$width <- 1
  par(mfrow = c(2, 4), mar = c(0, 0, 0, 0) + 0.5)
  for (cm in unique(node$cm)) {
    sub_graph <- subgraph(graph, vids = intersect(V(graph)$name,
                                                  node$subCluster[node$cm == cm]))
    sub_graph <- delete.vertices(sub_graph, V(sub_graph)[igraph::degree(sub_graph) ==
                                                           0])
    plot.igraph(sub_graph, layout = Layout, xlim = c(-1.2,
                                                     1.2), ylim = c(-1.2, 1.2), vertex.size = 50, vertex.label.cex = 5/8,
                vertex.label.color = "black", edge.curved = FALSE)
    title(cm, cex.main = 7/8, line = -0.5)
  }
}

### postprocessing
NMF_K6=readRDS("./AnaRes20250723/NMF_K4.rds")
NMF_K6_cor=readRDS("./AnaRes20250723/NMF_K4_cor_pair.rds")
#Mobj_keep@meta.data[,grep("^CMT",colnames(Mobj_keep@meta.data))]=NULL
#Mobj_keep@meta.data[,grep("^ScaleCMT",colnames(Mobj_keep@meta.data))]=NULL
NMF_K6_cor=NMF_K6_cor[abs(NMF_K6_cor$correlation)>=0.1 & NMF_K6_cor$pval<=0.25,]
#NMF_K6_cor[grepl("CNTN4",NMF_K6_cor$subCluster1) | grepl("CNTN4",NMF_K6_cor$subCluster2),]
CT_weight=NMF_K6@fit@W
colnames(CT_weight)=paste0("CMT0",1:4)
CMT_DEG=data.frame()
for(x in colnames(CT_weight)){
  tmp_score=CT_weight[,x]
  tmp_score=tmp_score[tmp_score>=0.02]
  tmp_vector=names(tmp_score)
  xxx=FindMarkers(Mobj_keep, ident.1 = tmp_vector, group.by = "subtypes015", only.pos = T)
  xxx=xxx[xxx$p_val_adj<=0.05 &
          xxx$avg_log2FC>=0.35 &
          xxx$pct.1-xxx$pct.2>=0.1,]
  xxx$gene=rownames(xxx)
  xxx$cluster=x
  CMT_DEG=rbind(CMT_DEG, xxx)
}

CMT_DEG=CMT_DEG[CMT_DEG$avg_log2FC>=1,]
for(x in unique(CMT_DEG$cluster)){
  Mobj_keep=AddModuleScore(Mobj_keep, features = list(CMT_DEG[CMT_DEG$cluster %in% x,]$gene),
                      ctrl = 500, name = paste0(x,".0"))
}
ModuleV=colnames(Mobj_keep@meta.data)[grepl("^CMT0", colnames(Mobj_keep@meta.data))]
for(m in ModuleV){
  mmm=paste0("Scale",m)
  Mobj_keep[[mmm]]=scale(Mobj_keep[[m]])[,1]
}
all_majorModule=c()
res.Module=Mobj_keep@meta.data[,grepl("^ScaleCMT0", colnames(Mobj_keep@meta.data))]
typy=colnames(res.Module) %>% gsub("Scale","",.)
for(i in 1:dim(res.Module)[1]){
  tidx=which(res.Module[i,] == max(res.Module[i,]), arr.ind = TRUE)[,2]
  all_majorModule=c(all_majorModule, typy[tidx])
}
Mobj_keep$CMT=gsub(".01$","",all_majorModule)

x="CMT04"
tmp_vector=names(CT_weight[,x][CT_weight[,x]>=0.02] %>% sort(., decreasing = TRUE))
df=NMF_K6_cor[NMF_K6_cor$subCluster1 %in% tmp_vector[1:10] &
           NMF_K6_cor$subCluster2 %in% tmp_vector[1:10],]
df=df[df$correlation>0,c("subCluster1","subCluster2","correlation")]
df_mtx=dcast(df, subCluster1 ~ subCluster2)
rownames(df_mtx)=df_mtx$subCluster1
df_mtx$subCluster1=NULL
df_mtx[is.na(df_mtx)]=0
df_mtx=as.matrix(df_mtx)
inmat=matrix_expand(df_mtx, identV = unique(c(rownames(df_mtx), colnames(df_mtx))))
#inmat[inmat<0.0000002]=0
inmat[inmat<0.1]=0
CairoPDF(paste0(OutputDir,"/Fig4.",x,"_top10_cor.pdf"), width = 10, height = 10)
CellChat::netVisual_circle(inmat, weight.scale = T, label.edge= T , vertex.size.max = 30, edge.curved = F, remove.isolate = F,
                           edge.weight.max = 0.5, edge.width.max = 10,  arrow.size = 0, alpha.edge = 0.35, margin = 0.02)
dev.off()
