# ============================================================================
# Plotting helpers (DoHeatmapSelf, fastpHeat, anno_pct, pick_group_color, splitDotPlot, plotpubr)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Source     : extracted from Analysis240525.r (L1311-1461, L1648-1661, L2297-2318, L2764-2793, L7216-7278, L8066-8097)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
DoHeatmapSelf=function (object, features = NULL, cells = NULL, group.by = "ident",
                      group.bar = TRUE, group.colors = NULL, disp.min = -2.5, disp.max = NULL,
                      slot = "scale.data", assay = NULL, label = TRUE, size = 5.5,
                      hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE,
                      lines.width = NULL, group.bar.height = 0.02, combine = FALSE) {

cells <- cells %||% colnames(x = object)
if (is.numeric(x = cells)) {
  cells <- colnames(x = object)[cells]
}
assay <- assay %||% DefaultAssay(object = object)
DefaultAssay(object = object) <- assay
features <- features %||% VariableFeatures(object = object)
features <- rev(x = unique(x = features))
disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                 yes = 2.5, no = 6)
possible.features <- rownames(x = GetAssayData(object = object,
                                               slot = slot))
if (any(!features %in% possible.features)) {
  bad.features <- features[!features %in% possible.features]
  features <- features[features %in% possible.features]
  if (length(x = features) == 0) {
    stop("No requested features found in the ", slot,
         " slot for the ", assay, " assay.")
  }
  warning("The following features were omitted as they were not found in the ",
          slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                           collapse = ", "))
}
in_data=as.matrix(x = t(x = GetAssayData(object = object,
                                         slot = slot)[features, cells, drop = FALSE]))
in_data=scale(in_data, center = TRUE, scale = TRUE)
data <- as.data.frame(x = in_data)
object <- suppressMessages(expr = StashIdent(object = object,
                                             save.name = "ident"))
group.by <- group.by %||% "ident"
groups.use <- object[[group.by]][cells, , drop = FALSE]
plots <- vector(mode = "list", length = ncol(x = groups.use))
for (i in 1:ncol(x = groups.use)) {
  data.group <- data
  group.use <- groups.use[, i, drop = TRUE]
  if (!is.factor(x = group.use)) {
    group.use <- factor(x = group.use)
  }
  names(x = group.use) <- cells
  if (draw.lines) {
    lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                              0.0025)
    placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                         lines.width), FUN = function(x) {
                                           return(RandomName(length = 20))
                                         })
    placeholder.groups <- rep(x = levels(x = group.use),
                              times = lines.width)
    group.levels <- levels(x = group.use)
    names(x = placeholder.groups) <- placeholder.cells
    group.use <- as.vector(x = group.use)
    names(x = group.use) <- cells
    group.use <- factor(x = c(group.use, placeholder.groups),
                        levels = group.levels)
    na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                            ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                         colnames(x = data.group)))
    data.group <- rbind(data.group, na.data.group)
  }
  lgroup <- length(levels(group.use))
  plot <- SingleRasterMap(data = data.group, raster = raster,
                          disp.min = disp.min, disp.max = disp.max, feature.order = features,
                          cell.order = names(x = sort(x = group.use)), group.by = group.use)
  if (group.bar) {
    default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
    if (!is.null(x = names(x = group.colors))) {
      cols <- unname(obj = group.colors[levels(x = group.use)])
    }
    else {
      cols <- group.colors[1:length(x = levels(x = group.use))] %||%
        default.colors
    }
    if (any(is.na(x = cols))) {
      cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
      cols <- Col2Hex(cols)
      col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                      start = 1, stop = 7)))))
      through <- length(x = default.colors)
      while (length(x = col.dups) > 0) {
        pal.max <- length(x = col.dups) + through
        cols.extra <- hue_pal()(pal.max)[(through +
                                            1):pal.max]
        cols[col.dups] <- cols.extra
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                        start = 1, stop = 7)))))
      }
    }
    group.use2 <- sort(x = group.use)
    if (draw.lines) {
      na.group <- RandomName(length = 20)
      levels(x = group.use2) <- c(levels(x = group.use2),
                                  na.group)
      group.use2[placeholder.cells] <- na.group
      cols <- c(cols, "#FFFFFF")
    }
    pbuild <- ggplot_build(plot = plot)
    names(x = cols) <- levels(x = group.use2)
    y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
    y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) +
      y.range * 0.015
    y.max <- y.pos + group.bar.height * y.range
    x.min <- min(pbuild$layout$panel_params[[1]]$x.range) +
      0.1
    x.max <- max(pbuild$layout$panel_params[[1]]$x.range) -
      0.1
    plot <- plot + annotation_raster(raster = t(x = cols[group.use2]),
                                     xmin = x.min, xmax = x.max, ymin = y.pos, ymax = y.max) +
      coord_cartesian(ylim = c(0, y.max), clip = "off") +
      scale_color_discrete(name = "Identity", na.translate = FALSE)
    if (label) {
      x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
      x.divs <- pbuild$layout$panel_params[[1]]$x.major %||%
        attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(),
             which = "pos")
      x <- data.frame(group = sort(x = group.use),
                      x = x.divs)
      label.x.pos <- tapply(X = x$x, INDEX = x$group,
                            FUN = function(y) {
                              if (isTRUE(x = draw.lines)) {
                                mean(x = y[-length(x = y)])
                              }
                              else {
                                mean(x = y)
                              }
                            })
      label.x.pos <- data.frame(group = names(x = label.x.pos),
                                label.x.pos)
      plot <- plot + geom_text(stat = "identity", data = label.x.pos,
                               aes_string(label = "group", x = "label.x.pos"),
                               y = y.max + y.max * 0.03 * 0.5, angle = angle,
                               hjust = hjust, size = size)
      plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                               y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) *
                                                                 size), clip = "off"))
    }
  }
  plot <- plot + theme(line = element_blank())
  plots[[i]] <- plot
}
if (combine) {
  plots <- wrap_plots(plots)
}
return(plots)
}

fastpHeat=function(zzz, rowcol="r"){
if(rowcol=="r"){
  zzzr=rowPercents(zzz)[,1:dim(zzz)[2]]
  pheatmap(zzzr, name = "pct.",
           cluster_rows = F, cluster_cols = F,
           display_numbers = zzzr)
}
if(rowcol=="c"){
  zzzc=colPercents(zzz)[1:dim(zzz)[1],]
  pheatmap(zzzc, name = "pct.",
           cluster_rows = F, cluster_cols = F,
           display_numbers = zzzc)
}
}

anno_pct = function(x, rc="row") {
  # x numbers matching row or columns
  max_x = max(x)
  text = paste0(sprintf("%.0f", x))
  cell_fun_pct = function(i) {
    pushViewport(viewport(xscale = c(0, max_x)))
    grid.roundrect(x = unit(1, "npc"), width = unit(x[i], "native"),
                   height = unit(1, "npc") - unit(4, "pt"),
                   just = "right", gp = gpar(fill = "#0000FF80", col = NA))
    grid.text(text[i], x = unit(1, "npc"), just = "right")
    popViewport()
  }
  AnnotationFunction(
    cell_fun = cell_fun_pct,
    var_import = list(max_x, x, text),
    which = rc,
    width = max_text_width(text)*1.25
  )
}

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

splitDotPlot=function(obj, group.by="cluster", groupV=NULL, split.by="group", splitV=NULL, featureV=NULL){

  library(scPioneer)

  if(!is.null(groupV)){
    obj=obj[, obj@meta.data[,group.by] %in% groupV]
    obj@meta.data[,group.by]=droplevels(obj@meta.data[,group.by])
  }
  if(!is.null(splitV)){
    obj=obj[, obj@meta.data[,split.by] %in% splitV]
    obj@meta.data[,split.by]=droplevels(obj@meta.data[,split.by])
  }
  groupV=levels(obj@meta.data[,group.by])
  splitV=levels(obj@meta.data[,split.by])

  targetGenes_select=featureV
  #targetGenes %>% group_by(cluster) %>% top_n(6, wt=avg_log2FC*(pct.1-pct.2)) %>% as.data.frame()
  #  xxx=scPioneer::DotPlot2(obj,
  #  xxx$data$group=rep(x, dim(xxx$data)[1])
  #}
  xxx=DotPlot(obj, features = unique(targetGenes_select), cols = hue_pal()(100), group.by = group.by, split.by = split.by)
  xxxdf=xxx$data
  xxxdf$group=gsub(paste0("^",paste(groupV, sep = "", collapse = "_|^"),"_"),
                   "", xxxdf$id)
  xxxdf$id=gsub(paste0("_",paste(splitV, sep = "", collapse = "$|_"),"$"),
                "", xxxdf$id)
  xxxdf$id=factor(xxxdf$id,levels = groupV)
  xxxdf$features.plot=factor(xxxdf$features.plot, levels = unique(targetGenes_select))
  xxxdf$group=factor(xxxdf$group, levels = splitV)

  col.min = -2.5; col.max = 2.5; scaleRun=TRUE; #dot.min = 0; dot.scale = 6;
avg.exp.scaled <- sapply(X = unique(x = xxxdf$features.plot),
                         FUN = function(x) {
                           data.use <- xxxdf[xxxdf$features.plot ==
                                               x, "avg.exp"]
                           if (scaleRun) {
                             data.use <- scale(x = log1p(data.use))
                             data.use <- MinMax(data = data.use, min = col.min,
                                                max = col.max)
                           }
                           else {
                             data.use <- log1p(x = data.use)
                           }
                           return(data.use)
                         })
avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
xxxdf$avg.exp.scaled <- avg.exp.scaled
  tplot=ggplot(xxxdf)+geom_point(aes(group, id, size=pct.exp, color=avg.exp.scaled))+
    facet_wrap(~features.plot, nrow = 1)+
    scale_color_gradientn(
      name = "Average\nExpression", colors = rev(scPioneer::scPalette_heatmap_bar())) +
    cowplot::theme_cowplot()+
    xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1), strip.text.x = element_text(angle = 75))

  return(tplot)
}

plotpubr=function(indf, gfacet=NULL, cmp=NULL, ...){

  splitdf1=ggpubr::compare_means(data = indf, formula = exp~group, group.by = gfacet, method = "t.test")
  stat1=ggpubr::stat_compare_means(data = indf, aes(label = ..p.signif..), vjust=0, bracket.size=0.5, hide.ns = FALSE, comparisons = cmp, method = "t.test")

  hhsplitbox1=ggpubr::ggviolin(indf, x = "group", y = "exp", width = 0.8, scale="width",
                               color="group", fill="group", trim = TRUE,
                               add = c("boxplot"),
                               add.params=list(shape=21, color=c("black"), fill="white", size=0.3, width=0.3),
                               short.panel.labs = TRUE,  #facet.by = gfacet, nrow=1,
                               ggtheme =theme_minimal())

  hhsplitbox1=hhsplitbox1+stat1

  return(hhsplitbox1)

}
