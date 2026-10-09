# ============================================================================
# Cellular co-occurrence modules with CoVarNet (frequency matrix, NMF factorisation, module plotting)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Identification of cellular modules"
# Source     : extracted from Analysis240525.r (L13268-13588, L13676-13790)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
########################## CoVarNet 202506 ###################

library(CoVarNet)
library(NMF)

# Preparing annotation files
# Filter out samples with cell numbers less than 100. Insufficient cell numbers can lead to errors in cell frequency.
meta <- Mobj_keep@meta.data #readRDS('rdata/merge_all_clean_meta.rds')
meta$SampleName=meta$sample
rt_sp <- names(table(meta$SampleName))[table(meta$SampleName) >= 100]
length(rt_sp)
meta <- meta[meta$SampleName %in% rt_sp, ]

# Cell subset frequency matrix

mat_fq_raw <- freq_calculate(meta,major_anno_col = "refine_ct03",
                             sub_anno_col = "subtypes015",
                             sampleID_col = "SampleName")
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize="minmax")
mat_fq_norm=mat_fq_norm[,colSums(mat_fq_norm)>0]
mat_fq_raw=mat_fq_raw[, c("refine_ct03",colnames(mat_fq_norm))]

#Performing NMF decomposition on frequency matrix
res <- nmf(as.matrix(mat_fq_norm), rank = 2:20, method = "nsNMF", seed = rep(123456, 6), .options = "vp")
plot(res)

## plot pipeling
plotting <- function(mat_fq_norm, mat_fq_raw, K, outdir){
  NMF_res <- nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
  saveRDS(NMF_res, paste0(outdir, '/NMF_K',K, '.rds'))
  #K=8

  pdf(paste0(outdir, '/NMF_K',K, '_weight_htmap.pdf'), 6,48)
  gr.weight_all(NMF_res,wid=5,hei=40)
  dev.off()

  pdf(paste0(outdir, '/NMF_K',K, '_weight_top20.pdf'), 20,30)
  p <- gr.weight_top(NMF_res,num=30)
  print(p)
  dev.off()

  pdf(paste0(outdir, '/NMF_K',K,'_distribution_htmap.pdf'))
  gr.distribution(NMF_res,meta=meta,group="Organ",sampleID_col='SampleName',extra_anno_col='Location_short')
  dev.off()

  ### digit: 3
  mat_fq_raw[,2:dim(mat_fq_raw)[2]]=round(mat_fq_raw[,2:dim(mat_fq_raw)[2]],digits = 3)
  cor_pair<-pair_correlation(
    mat_fq_raw[highly_variable_samples,],
    method="spearman")
  saveRDS(cor_pair, paste0(outdir, '/NMF_K',K,'_cor_pair.rds'))

  top_n=8
  network<-cm_network(NMF_res,cor_pair,top_n,corr=0.05,fdr=0.25)
  each<-network$each
  pdf(paste0(outdir, '/NMF_K',K,'_network.pdf'),10,5)
  Xgr.igraph_each(each,Layout=layout_in_circle)
  dev.off()
}

### running
dir.create("./AnaRes20250701/")
outdir <- './AnaRes20250701/'
plotting(mat_fq_norm, mat_fq_raw, 3, outdir)
plotting(mat_fq_norm, mat_fq_raw, 4, outdir)
plotting(mat_fq_norm, mat_fq_raw, 6, outdir)
plotting(mat_fq_norm, mat_fq_raw, 8, outdir)
plotting(mat_fq_norm, mat_fq_raw, 12, outdir)

### refine 0723
dir.create("./AnaRes20250723/")
outdir <- './AnaRes20250723/'
plotting(mat_fq_norm[rownames(mat_fq_norm) %!in% c("Prolif."),],
         mat_fq_raw[rownames(mat_fq_norm) %!in% c("Prolif."),], 4, outdir)

freq_calculate <- function(meta, major_anno_col, sub_anno_col, sampleID_col) {
  # 确保列名存在
  if (!all(c(major_anno_col, sub_anno_col, sampleID_col) %in% colnames(meta))) {
    stop("One or more specified columns do not exist in the metadata")
  }

  # 获取唯一的subCluster对应的majorCluster信息
  meta_sc <- meta[!duplicated(meta[[sub_anno_col]]),
                  c(major_anno_col, sub_anno_col)]

  # 设置行名并排序
  rownames(meta_sc) <- meta_sc[[sub_anno_col]]
  meta_sc <- meta_sc[order(meta_sc[[sub_anno_col]]), ]

  # 创建计数矩阵
  mat_ct <- table(meta[[sub_anno_col]], meta[[sampleID_col]])
  mat_ct <- matrix(data = mat_ct, ncol = ncol(mat_ct),
                   dimnames = dimnames(mat_ct)) %>%
    as.data.frame()

  # 匹配顺序
  id <- match(rownames(mat_ct), meta_sc[[sub_anno_col]])
  mat_ct <- data.frame(majorCluster = factor(meta_sc[[major_anno_col]][id]),
                       mat_ct, check.names = FALSE)

  # 计算频率
  mat_fq <- mat_ct %>%
    group_by(majorCluster) %>%
    mutate_if(is.numeric, .funs = list(~./sum(.))) %>%
    as.data.frame()

  # 处理NA值
  mat_fq[is.na(mat_fq)] <- 0
  rownames(mat_fq) <- rownames(mat_ct)

  # 重命名majorCluster列以反映实际使用的列名
  colnames(mat_fq)[1] <- major_anno_col

  return(mat_fq)
}

gr.weight_all <- function(nmf_res,wid=5,hei=15){
  w <- basis(nmf_res)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
  }
  ht1 <- ComplexHeatmap::Heatmap(matrix = w, name = "weight",
                                 circlize::colorRamp2(seq(0, 0.1, length.out = 5), viridis::viridis(n = 5)),
                                 width = unit(wid, "cm"), height = unit(hei, "cm"), row_names_gp = grid::gpar(fontsize = 4.5),
                                 column_names_gp = grid::gpar(fontsize = 8), clustering_method_rows = "ward.D2")
  ComplexHeatmap::draw(ht1)
}

gr.distribution <- function(nmf_res, sampleID_col, extra_anno_col = NULL, ...) {
  var_args <- list(...)
  meta <- if (!is.null(var_args[["meta"]])) var_args[["meta"]] else NULL
  group <- if (!is.null(var_args[["group"]])) var_args[["group"]] else NULL

  h <- scoef(nmf_res)

  # 设置默认行名（成分名）
  if (is.null(rownames(h))) {
    rownames(h) <- sprintf("CM%02d", 1:nrow(h))
  }

  # 排序成分
  h <- h[order(rownames(h)), ]

  # 获取每个样本的主要成分
  id <- apply(h, MARGIN = 2, FUN = which.max)

  # 处理meta数据
  if (is.null(meta)) {
    meta_sp <- data.frame(
      CMT = sprintf("CMT%02d", 1:nrow(h))[id],
      sampleID = colnames(h),
      stringsAsFactors = FALSE
    )
    rownames(meta_sp) <- meta_sp$sampleID
    meta_sp <- meta_sp[order(meta_sp$CMT), ]
  } else {
    meta_sp <- meta[!duplicated(meta[, sampleID_col]), ]
    rownames(meta_sp) <- meta_sp[, sampleID_col]
    meta_sp <- meta_sp[colnames(h), ]
    meta_sp$CMT <- sprintf("CMT%02d", 1:nrow(h))[id]

    # 排序
    if (is.null(group)) {
      meta_sp <- meta_sp[order(meta_sp$CMT), ]
    } else {
      meta_sp <- meta_sp[order(meta_sp$CMT, meta_sp[, group]), ]
    }
  }

  # 确保热图数据与meta顺序一致
  h <- h[, rownames(meta_sp)]

  # 设置CMT颜色（固定使用Set3调色板）
  CMTnames <- gsub("CM", "CMT", rownames(h))
  cmt_colors <- if (nrow(h) <= 12) {
    RColorBrewer::brewer.pal(nrow(h), "Set3")
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(h))
  }
  cmt_colors <- cmt_colors[1:nrow(h)]
  names(cmt_colors) <- CMTnames

  # 构建注释列表和颜色列表
  anno_list <- list(CMT = meta_sp$CMT)
  col_list <- list(CMT = cmt_colors)

  # 为group列分配颜色（使用Paired调色板）
  if (!is.null(group) && group %in% colnames(meta_sp)) {
    group_vals <- meta_sp[, group]
    if (is.factor(group_vals) || is.character(group_vals)) {
      unique_groups <- unique(group_vals)
      n_groups <- length(unique_groups)

      # 为group选择专用调色板（Paired）
      group_colors <- if (n_groups <= 12) {
        RColorBrewer::brewer.pal(max(3, n_groups), "Paired")[1:n_groups]
      } else {
        colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_groups)
      }
      names(group_colors) <- unique_groups

      anno_list[["Group"]] <- group_vals
      col_list[["Group"]] <- group_colors
    }
  }

  # 添加额外注释列（使用Accent调色板）
  if (!is.null(extra_anno_col) && extra_anno_col %in% colnames(meta_sp)) {
    extra_vals <- meta_sp[, extra_anno_col]
    if (is.factor(extra_vals) || is.character(extra_vals)) {
      unique_extra <- unique(extra_vals)
      n_extra <- length(unique_extra)

      # 为额外列选择专用调色板（Accent）
      extra_colors <- if (n_extra <= 8) {
        RColorBrewer::brewer.pal(max(3, n_extra), "Set2")[1:n_extra]
      } else {
        colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(n_extra)
      }
      names(extra_colors) <- unique_extra

      anno_list[[extra_anno_col]] <- extra_vals
      col_list[[extra_anno_col]] <- extra_colors
    }
  }

  # 创建热图注释
  ann_top <- ComplexHeatmap::HeatmapAnnotation(
    df = anno_list,
    col = col_list,
    border = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 6),
    simple_anno_size = unit(3, "mm") ,
    annotation_legend_param = list(
      legend_direction = "horizontal",
      ncol = 2,
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6),
      grid_width = unit(3, "mm"),
      grid_height = unit(2, "mm")
    )
  )

  # 创建热图（其余部分保持不变）
  ht <- ComplexHeatmap::Heatmap(
    h,
    width = unit(3, "in"),
    height = unit(1.5, "in"),
    name = "CM abundance",
    border = TRUE,
    show_column_names = FALSE,
    row_names_side = "left",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = "CM abundance in samples",
    column_title_gp = grid::gpar(fontsize = 8),
    col = circlize::colorRamp2(seq(0, 0.5, length.out = 4), viridis::viridis(4)),
    top_annotation = ann_top,
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      color_bar = "continuous",
      title_position = "lefttop",
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6),
      grid_height = unit(3, "mm")
    )
  )

  # 绘制热图（其余部分保持不变）
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )

  # 添加装饰线（其余部分保持不变）
  dup <- (which(!duplicated(meta_sp$CMT)) - 1)
  fract <- dup/nrow(meta_sp)
  width <- c(fract[-1], 1) - fract
  ComplexHeatmap::decorate_heatmap_body("CM abundance", {
    grid::grid.rect(
      x = unit(fract, "native"),
      y = unit(1 - (0:(length(dup) - 1))/length(dup), "native"),
      width = unit(width, "native"),
      height = unit(1/length(unique(meta_sp$CMT)), "native"),
      hjust = 0,
      vjust = 1,
      gp = grid::gpar(col = "white", fill = NA, lty = 1, lwd = 2)
    )
  })
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
