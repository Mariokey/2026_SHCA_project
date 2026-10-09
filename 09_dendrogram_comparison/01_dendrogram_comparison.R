# ============================================================================
# Organ developmental dendrogram comparison and Robinson-Foulds distances
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Dendrogram comparison"
# Source     : extracted from tree_compare.r (L1-486)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
# Load required packages
if (!requireNamespace("ggtree", quietly = TRUE)) {
  stop("ggtree package is required: install.packages('ggtree')")
}
if (!requireNamespace("phangorn", quietly = TRUE)) {
  stop("phangorn package is required: install.packages('phangorn')")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr package is required: install.packages('dplyr')")
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("需要ape包，请先安装：install.packages('ape')")
}
library(ggtree)
library(phangorn)
library(dplyr)
library(ape)

### 20250812

reftree_df=read.xlsx("./organ_combinations_with_relationship_updated.xlsx")
reftree_df$organ1=gsub("Gallbladder","GBD",reftree_df$organ1)
reftree_df$organ2=gsub("Gallbladder","GBD",reftree_df$organ2)
reftree_df$organ1=gsub("Ovary","OV",reftree_df$organ1)
reftree_df$organ2=gsub("Ovary","OV",reftree_df$organ2)
reftree_df$organ1=gsub("\\s+","",reftree_df$organ1)
reftree_df$organ2=gsub("\\s+","",reftree_df$organ2)
reftree=tidyr::pivot_wider(reftree_df, id_cols = "organ1", names_from  = "organ2", values_from = "relationship")
reftree_mat=reftree[,2:dim(reftree)[2]] %>% as.matrix()
rownames(reftree_mat)=reftree[[1]] %>% as.vector()
#xxx=pheatmap(reftree_mat, name = "ref_dis", clustering_method = "wald.D")

allct_color=hue_pal()(length(rownames(reftree_mat)))
names(allct_color)=rownames(reftree_mat)

cmptree_list=list()
for(x in 1:30){

  #x=43
  x_ct=levels(Mobj_keep$subtypes015)[x+31]
  tmpdata=subset(Mobj_keep, subtypes015 %in% x_ct)

xxx_tbl=table(tmpdata$Location_short)/dim(tmpdata)[2]
print(paste0(x_ct,": ",paste(names(xxx_tbl[xxx_tbl>=0.05]), collapse=" ")))
if(x_ct %in% c("Hepatocytes")){
  tmpdata=subset(tmpdata, Location_short %in% names(xxx_tbl[xxx_tbl>=0.1]))
}else if(x_ct %in% c("Tuft_Cells","CNTN4_NEL_Cells")){
  tmpdata=subset(tmpdata, Location_short %in% names(xxx_tbl[xxx_tbl>=0.005]))
}else{
  tmpdata=subset(tmpdata, Location_short %in% names(xxx_tbl[xxx_tbl>=0.05]))
}
tmpdata$Location_short=droplevels(tmpdata$Location_short)

if(length(levels(tmpdata$Location_short))>=3){
tmpdata=FindVariableFeatures(tmpdata, verbose = 0)
LS_Gene <- VariableFeatures(object = tmpdata)
tdata <- subset(tmpdata,features = LS_Gene)
Idents(tdata)="Location_short"
tdata <- AverageExpression(tdata,
                           return.seurat = TRUE)
Mat_Tmp <- as.matrix(tdata@assays$RNA$data)
hc= hclust(dist(t(Mat_Tmp)))
ct_color=allct_color[colnames(Mat_Tmp)]

cmptree_list[[x_ct]]=compare_trees(reftree_mat, hc, ct_color,
                               method = c("topology"))
}

}

pdf("./tree_compare.pdf")
tmpdf_name=c()
tmpdf_sim=c()
for(x in 1:23){
  tmpdf_name=c(tmpdf_name, names(cmptree_list)[x])
  tmpdf_sim=c(tmpdf_sim, round(cmptree_list[[x]]$similarity, digits = 3))
  title=paste0(names(cmptree_list)[x],"_sim: ",
               round(cmptree_list[[x]]$similarity, digits = 3))
  grid::grid.newpage()
  grid.text(title,x = (0.5), y = (0.6))
  }
dev.off()
tmpdf=data.frame(celltype=tmpdf_name, similarity=tmpdf_sim)

### 202509
cmptree_CTlist=list()
for(x in c(1, 3, 4, 6, 7, 26, 27)){   ### 1 3 4 6 7 26 27

  #x=1
  x_ct=levels(Mobj_keep$refine_ct03)[x]
  tmpdata=subset(Mobj_keep, refine_ct03 %in% x_ct)

  xxx_tbl=table(tmpdata$Location_short)
  print(paste0(x_ct,": ",paste(names(xxx_tbl[xxx_tbl>=10]), collapse=" ")))
  #}
  tmpdata=subset(tmpdata, Location_short %in% names(xxx_tbl[xxx_tbl>=10]))
  tmpdata$Location_short=droplevels(tmpdata$Location_short)

  if(length(levels(tmpdata$Location_short))>=3){
    tmpdata=FindVariableFeatures(tmpdata, verbose = 0)
    LS_Gene <- VariableFeatures(object = tmpdata)
    tdata <- subset(tmpdata,features = LS_Gene)
    Idents(tdata)="Location_short"
    tdata <- AverageExpression(tdata,
                               return.seurat = TRUE)
    Mat_Tmp <- as.matrix(tdata@assays$RNA$data)
    hc= hclust(dist(t(Mat_Tmp)))
    ct_color=allct_color[colnames(Mat_Tmp)]

    cmptree_CTlist[[x_ct]]=compare_trees(reftree_mat, hc, ct_color,
                                       method = c("topology"))
  }

}

for(x in names(cmptree_CTlist)){
  runFig(paste0(FinalFig,"/", x, "_treeplot"), cmptree_CTlist[[x]]$tree_plot, 9, 9)
}

xdf=data.frame(celltypes=c("B","T","NK","Mono","Mph","EC","FB","Tuft_Cells","CNTN4_NEL_Cells","Endocrine_Cells","Colonocytes","Goblet_Cells","Prezymogenic_Cells","Melanocytes","Multiciliated_Cells","Papillary_Epithelial_Cells"),
                devIndex=c(0.086,0.122,0.083,0.121,0.162,0.189,0.075,0.267,0.208,0.333,1,0.25,1,0.5,1,1))
lollipop_plot(xdf, "celltypes","devIndex", dot_color = allct_color)

### test
phylo2 <- ape::as.phylo(hc)

reftree_mat_sub=reftree_mat
overlap_lab=intersect(hc$labels,colnames(reftree_mat_sub))
reftree_mat_sub=reftree_mat_sub[overlap_lab,overlap_lab]
refTree_hc=hclust(RFLPtools::sim2dist(reftree_mat_sub, maxSim = 6))
phylo1 <- ape::as.phylo(refTree_hc) #tree1$phylo

splits1 <- ape::prop.part(phylo1)
splits2 <- ape::prop.part(phylo2)
unique_splits1 <- setdiff(seq_along(splits1), phangorn::matchSplits(splits1, splits2))
unique_splits2 <- setdiff(seq_along(splits2), phangorn::matchSplits(splits1, splits2))
split_sizes1 <- sapply(splits1[unique_splits1], length)
split_sizes2 <- sapply(splits2[unique_splits2], length)

all_splits <- c(
  setNames(split_sizes1, paste0("TreeRef_split_", unique_splits1)),
  setNames(split_sizes2, paste0("TreeQry_split_", unique_splits2))
)
result$different_clades <- sort(all_splits, decreasing = TRUE)
result$most_influential_clades <- names(result$different_clades)[1:min(6, length(result$different_clades))]

### tree compare
compare_trees <- function(refTree_mat, inTree_hc, ct_color, method = c("topology")) {

  ### reference subtree
  print(setdiff(inTree_hc$labels,colnames(refTree_mat)))
  overlap_lab=intersect(inTree_hc$labels,colnames(refTree_mat))
  reftree_mat_sub=refTree_mat[overlap_lab,overlap_lab]
  refTree_hc=hclust(RFLPtools::sim2dist(reftree_mat_sub, maxSim = 6))
  phylo1 <- ape::as.phylo(refTree_hc) #tree1$phylo
  ### input subtree
  phylo2 <- ape::as.phylo(inTree_hc) #tree2$phylo

  # 检查树是否具有相同的物种集合
  if (!setequal(phylo1$tip.label, phylo2$tip.label)) {
    warning("两个树的末端标签不完全相同，比较结果可能不准确")
  }

  # 标准化树的结构，确保它们可以比较
  phylo1 <- ape::multi2di(phylo1)  # 处理多分叉节点
  phylo2 <- ape::multi2di(phylo2)

  method <- match.arg(method)
  result <- list()

  # 计算整体相似性指标并识别影响最大的节点/分支
  if (method %in% c("rf", "wrf")) {
    # Robinson-Foulds距离比较
    rf_args <- list(phylo1, phylo2, normalize = TRUE, details = TRUE)
    if (method == "wrf") {
      rf_args$weighted <- TRUE
    }

    rf_result <- do.call(phangorn::RF.dist, rf_args)
    result$distance <- rf_result$distance
    result$similarity <- 1 - rf_result$distance
    result$method <- ifelse(method == "rf", "Robinson-Foulds", "Weighted Robinson-Foulds")
    result$different_clades <- rf_result$clades

    # 找到影响最大的节点
    if (!is.null(phylo1$node.label) && !is.null(phylo2$node.label)) {
      node_mapping <- phangorn::matchNodes(phylo1, phylo2)
      valid_matches <- !is.na(node_mapping[, 2])

      if (any(valid_matches)) {
        # 提取匹配节点的支持度并计算差异
        support1 <- as.numeric(phylo1$node.label[node_mapping[valid_matches, 1]])
        support2 <- as.numeric(phylo2$node.label[node_mapping[valid_matches, 2]])
        support_diff <- abs(support1 - support2)
        names(support_diff) <- paste0("Node_", node_mapping[valid_matches, 1], "_vs_", node_mapping[valid_matches, 2])

        # 按差异排序并获取影响最大的节点
        result$node_support_differences <- sort(support_diff, decreasing = TRUE)
        result$most_influential_nodes <- names(result$node_support_differences)[1:min(5, length(result$node_support_differences))]
      } else {
        warning("无法找到两个树之间的匹配节点进行支持度比较")
      }
    } else {
      warning("一个或两个树对象缺少节点标签，无法计算节点支持度差异")
    }
  } else if (method == "branch_length") {
    # 分支长度比较
    edge_mapping <- phangorn::matchEdges(phylo1, phylo2)
    valid_matches <- !is.na(edge_mapping[, 2])

    if (any(valid_matches)) {
      # 提取匹配分支的长度并计算差异
      len1 <- phylo1$edge.length[edge_mapping[valid_matches, 1]]
      len2 <- phylo2$edge.length[edge_mapping[valid_matches, 2]]
      len_diff <- abs(len1 - len2)
      rel_diff <- len_diff / pmax(len1, len2)  # 相对差异

      # 存储结果
      result$mean_absolute_diff <- mean(len_diff, na.rm = TRUE)
      result$mean_relative_diff <- mean(rel_diff, na.rm = TRUE)
      result$similarity <- 1 - result$mean_relative_diff
      result$method <- "Branch Length Comparison"

      # 按相对差异排序并获取影响最大的分支
      names(rel_diff) <- paste0("Edge_", edge_mapping[valid_matches, 1], "_vs_", edge_mapping[valid_matches, 2])
      result$branch_length_differences <- sort(rel_diff, decreasing = TRUE)
      result$most_influential_branches <- names(result$branch_length_differences)[1:min(5, length(result$branch_length_differences))]
    } else {
      warning("无法找到两个树之间的匹配分支进行长度比较")
    }
  } else if (method == "topology") {
    # 拓扑结构比较（忽略分支长度）
    topo1 <- phylo1
    topo1$edge.length <- NULL
    topo2 <- phylo2
    topo2$edge.length <- NULL

    # 检查拓扑结构是否完全相同
    result$topology_identical <- ape::all.equal.phylo(topo1, topo2) == TRUE
    result$method <- "Topology Comparison"

    if (!result$topology_identical) {
      # 使用treedist获取拓扑差异
      td_result <- phangorn::RF.dist(topo1, topo2, normalize = TRUE, check.labels = TRUE)
      result$distance <- td_result
      result$similarity <- 1 - result$distance / max(td_result, 1)

      # 识别same分支
      splits1 <- ape::prop.part(topo1)
      splits2 <- ape::prop.part(topo2)

      unique_splits1 <- setdiff(seq_along(splits1), phangorn::matchSplits(splits1, splits2))
      unique_splits2 <- setdiff(seq_along(splits2), phangorn::matchSplits(splits2, splits1))
      split_sizes1 <- sapply(splits1[unique_splits1], length)
      split_sizes2 <- sapply(splits2[unique_splits2], length)

      all_splits <- c(
        setNames(split_sizes1, paste0("TreeRef_split_", unique_splits1)),
        setNames(split_sizes2, paste0("TreeQry_split_", unique_splits2))
      )
      result$different_clades <- sort(all_splits, decreasing = TRUE)
      result$most_influential_diffClades <- names(result$different_clades)[1:min(6, length(result$different_clades))]

      unique_splits1 <- intersect(seq_along(splits1), phangorn::matchSplits(splits1, splits2))
      unique_splits2 <- intersect(seq_along(splits2), phangorn::matchSplits(splits2, splits1))
      if(length(unique_splits1)>0 & length(unique_splits2)>0){
      split_sizes1 <- sapply(splits1[unique_splits1], length)
      split_sizes2 <- sapply(splits2[unique_splits2], length)

      all_splits <- c(
        setNames(split_sizes1, paste0("TreeRef_split_", unique_splits1)),
        setNames(split_sizes2, paste0("TreeQry_split_", unique_splits2))
      )
      result$same_clades <- sort(all_splits, decreasing = TRUE)
      result$most_influential_sameClades <- names(result$same_clades)[1:min(6, length(result$same_clades))]
      }
    } else {
      result$distance <- 0
      result$similarity <- 1
    }
  }

  # 添加基本信息
  result$n_tips_ref <- phylo1$tip.label
  result$n_tips_qry <- phylo2$tip.label
  result$n_nodes_ref <- phylo1$Nnode
  result$n_nodes_qry <- phylo2$Nnode

  result$tree_plot=drawTreeV0(phylo1, ct_color, layout_hv="v")+
                   drawTreeV0(phylo2, ct_color, layout_hv="v")

  highlight_influential_clades(phylo2, comparison_result = cmptree_CTlist$B,
                               type="diff", tree_id = "TreeRef")

  return(result)
}

highlight_influential_clades <-
  function(phylo_tree, type="diff", comparison_result,
           tree_id = "TreeRef") {

    phylo_tree=phylo2
    type="diff"
    comparison_result= cmptree_CTlist$B
    tree_id = "TreeQry"

  # Extract the splits (clades) from the original tree
  splits <- ape::prop.part(phylo_tree)

  if(type=="diff"){
  # Get influential clade names for the specified tree (e.g., "Tree1_split_3")
  influential_clades <- comparison_result$most_influential_diffClades
  target_clades <- grep(paste0("^", tree_id, "_split_"), influential_clades, value = TRUE)
  }else{
    influential_clades <- comparison_result$most_influential_sameClades
    target_clades <- grep(paste0("^", tree_id, "_split_"), influential_clades, value = TRUE)
  }

  if (length(target_clades) == 0) {
    warning(paste("No influential clades found for", tree_id))
    return(ggtree::ggtree(phylo_tree))
  }

  # Extract split IDs (e.g., 3 from "Tree1_split_3")
  split_ids <- as.integer(sub(paste0(tree_id, "_split_"), "", target_clades))

  # Get tip labels for each influential clade
  clade_tips <- lapply(split_ids, function(id) {
    splits[[id]]  # Returns tip labels in this clade
  })
  names(clade_tips) <- target_clades

  # Create a ggtree plot
  p <- ggtree::ggtree(phylo_tree)

  # Highlight each influential clade with a different color
  colors <- scales::hue_pal()(length(clade_tips))
  names(colors) <- target_clades

  tip_add=max(phylo_tree$edge.length)/50   ###0.5
  label_add=tip_add*3

  for (i in seq_along(clade_tips)) {
    i=1

    clade_name <- names(clade_tips)[i]
    tips <- clade_tips[[i]]
    color <- colors[i]

    # Highlight the clade with a background rectangle
    p <- p +
      ggtree::geom_cladelabel(
      node = phytools::findMRCA(phylo_tree, tips),  # Find most recent common ancestor node
      label = "", #clade_name,
      color = color,
      offset = label_add*4,
      offset.text = 1,
      angle=90,
      hjust=0.5,
      barsize = 1.5,
      fontsize = 3.5
    ) +
      # Color the tips in the clade
      ggtree::geom_tippoint(
        data = subset(p$data, label %in% tips),
        color = color,
        size = 3
      ) +
      geom_tiplab(angle=0, hjust=0, offset=label_add/2, show.legend=FALSE) +
      theme_tree(plot.margin=margin(10,100,10,10))+
      theme(legend.position="none")
  }

  return(p)
}

drawTreeV0=function(phylo_tree, ct_color, layout_hv="h", addtree=NULL, ...){

### NOTE: phylo_tree can be phylo or hc

p <- ggtree(phylo_tree, linetype='solid')
d <- data.frame(label = phylo_tree$tip.label,
                cyl = phylo_tree$tip.label)

tip_add=max(phylo_tree$edge.length)/50   ###0.5
label_add=tip_add*3     ### 1.5

if(!is.null(addtree)){

  phylo11=phylo_tree
  phylo11$edge.length=NULL
  phylo22=addtree
  phylo22$edge.length=NULL
  p1 <- ggtree(phylo11)
  p2 <- ggtree(phylo22)

  dd1 <- p1$data
  dd2 <- p2$data

  ## reverse x-axis and
  ## set offset to make the tree in the right hand side of the first tree
  dd2$x <- max(dd2$x) - dd2$x + max(dd1$x) + 10

  pp <- p1 + geom_tiplab() + geom_tree(data=dd2) + geom_tiplab(data = dd2, hjust=1)

  dd <- bind_rows(dd1, dd2) %>% filter(!is.na(label))

  p=pp + geom_line(aes(x, y, group=label), data=dd, color='lightgrey') +
    geom_treescale() + geom_treescale(x=8)

}else{
if(layout_hv=="h"){
p = p %<+% d +
  layout_dendrogram() +
  geom_tippoint(aes(fill=factor(cyl), x=x+tip_add),
                size=8, shape=21, color='white') +
  geom_tiplab(angle=45, hjust=1, offset=-1*(2*label_add), show.legend=FALSE) +
  scale_fill_manual(values = ct_color)+
  theme_tree(plot.margin=margin(50,50,80,50))+
  theme(legend.position="none")
}else{
p = p %<+% d +
  layout_rectangular() +
  geom_tippoint(aes(fill=factor(cyl), x=x+tip_add),
                size=4, shape=21, color='white') +
  geom_tiplab(angle=0, hjust=0, offset=label_add, show.legend=FALSE) +
  scale_fill_manual(values = ct_color)+
  theme_tree(plot.margin=margin(50,50,80,50))+
  theme(legend.position="none")
}
}
return(p)

}
