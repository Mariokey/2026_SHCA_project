# ============================================================================
# System-wise epithelial subtype refinement (subsysList)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L10219-10540)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
#################################################################################
x_abV=c("x01","x02","x03","x04")

plot_tmplist=list()
for(xi in 1:length(names(subsysList))){

  x=names(subsysList)[xi]

  x_a=eval(parse(text=paste0(x_abV[xi],"v")))
  x_b=eval(parse(text=paste0(x_abV[xi],"v_name")))

  Idents(subsysList[[x]])="refine_subtypes01"
  tmpdata=DotPlot(subsysList[[x]], features = names(x_a), group.by = "refine_subtypes01")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
  tmpdata$data$Descrip=x_b[as.vector(x_a[tmpdata$data$features.plot])] %>% as.vector()
  tmpdata$data$Descrip=factor(tmpdata$data$Descrip, levels = unique(as.vector(x_b)))
  tmpxy=tmpdata$mapping
  tmpxy$x=expr(`Descrip`)
  tmpdata$mapping=tmpxy

  plot_tmplist[[x]]=tmpdata
}

Sys2CT=rbind(unique(subsysList$RS@meta.data[,c("System","refine_subtypes01")]),
             unique(subsysList$US@meta.data[,c("System","refine_subtypes01")]),
             unique(subsysList$DS@meta.data[,c("System","refine_subtypes01")]),
             unique(subsysList$SKIN@meta.data[,c("System","refine_subtypes01")]))
rownames(Sys2CT)=NULL
colnames(Sys2CT)=c("Sys","CT")
Sys2CT=Sys2CT[Sys2CT$CT %!in% as.vector(Sys2CT$CT[duplicated(Sys2CT$CT)]),]
Sys2CT_v=Sys2CT$Sys
names(Sys2CT_v)=Sys2CT$CT

GO_score_df=rbind(plot_tmplist$RS$data %>% dplyr::mutate(Sys="RS"),
       plot_tmplist$US$data %>% dplyr::mutate(Sys="US"),
       plot_tmplist$DS$data %>% dplyr::mutate(Sys="DS"),
       plot_tmplist$SKIN$data %>% dplyr::mutate(Sys="SKIN"))
GO_score_df$SysCT=paste0(GO_score_df$Sys,".",GO_score_df$id)
GO_score_df=GO_score_df[!duplicated(paste0(GO_score_df$Descrip,GO_score_df$SysCT)),]

GO_score_mtx=dcast(GO_score_df[,c("Descrip","SysCT","avg.exp")],formula = "Descrip~SysCT") # %>% pivot_wider(id_cols = "Descrip", names_from = "SysCT", values_from = "avg.exp")
GO_score_mtx00=GO_score_mtx[,-1] %>% as.matrix()
rownames(GO_score_mtx00)=GO_score_mtx$Descrip #GO_score_mtx$Description
#GO_score_mtx00[is.na(GO_score_mtx00)]=min(GO_score_df$AvgScore)

pheatmap(
         GO_score_mtx00[rownames(GO_score_mtx00) %in% GO_DATA$PATHID2NAME[c("GO:0016055", "GO:0001701", "GO:0019827", "GO:0048863", "GO:0030154", "GO:0048699", "GO:0048666", "GO:0008033", "GO:0030425", "GO:0007155", "GO:0007411")] %>% as.vector(),],
         scale = F,
         cluster_rows = T,
         cluster_cols = T,
         color = viridis_pal()(1000))

##############################
#scPioneer::DimPlot_idx
DimPlot(#subsysListNEW$RS,
        subsysListNEW$US,
        cols=berryFunctions::addAlpha(mycolortmp[levels(subsysListNEW$US$refine_subtypes01)],0.5),
        group.by = "refine_subtypes01", label = T, repel = T)+ #, idx.sep = "  ", prefix.index = NULL)+
  NoLegend()+NoTitle+
  NoAxes()

#xxx=read_tsv("./AnaRes20250110/Location_short_subtyperove_main_age.txt")
Imm2Loc_roe_long$name=factor(Imm2Loc_roe_long$name, levels = colnames(xxx)[-1])
yyy=ggplot(Imm2Loc_roe_long,aes(name,log(1+value),fill=class))+geom_violin(scale = "width", trim = FALSE)+
  theme_cowplot()+
  xlab("")+ylab("log(ro/e+1)")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
yyy=yyy+stat_summary(fun.y = median, geom = "point", shape=95, size = 3, color = "darkred", position = position_dodge(width = 0.9))

#xxx=read.xlsx("./SysCT_DescripGO_MeanScore.xlsx")
#xxx$avg.exp=xxx$AvgScore
xxx=read.xlsx("./SysCT_DescripGO_MeanScore_all156.xlsx")
ct2go_json=convert_to_cytoscape_json(xxx)
visualize_network_with_weights(xxx[xxx$avg.exp>=median(xxx$avg.exp),], threshold = 0.02)

#xxx=read.xlsx("./SysCT_DescripGO_MeanScore_all156.xlsx")

visualize_network_with_weights <- function(df,threshold=0.2) {

  library(igraph)
  library(ggraph)
  library(dplyr)

  # Create a data frame of edges with weights
  edges <- data.frame(
    from = df$SysCT,
    to = df$Descrip,
    weight = df$avg.exp  # Use AvgScore as edge weight
  )

  # Filter edges based on the threshold
  edges <- edges %>% filter(weight >= threshold)

  # Create an igraph object
  g <- graph_from_data_frame(edges, directed = TRUE)

  # Set edge weights
  E(g)$weight <- edges$weight

  # Create a data frame for node attributes
  nodes <- data.frame(name = V(g)$name)

  # Add a column to distinguish SysCT and Description
  nodes <- nodes %>%
    mutate(type = ifelse(grepl("^RS\\.", name), "RS",
                         ifelse(grepl("^DS\\.", name), "DS",
                                ifelse(grepl("^US\\.", name), "US",
                                      ifelse(grepl("^SKIN\\.", name), "SKIN", "Description_GO")))))

  # Add node attributes to the igraph object
  V(g)$type <- nodes$type

  # Plot the network
  ggraph(g, layout = "fr") +  # Use Fruchterman-Reingold layout
    geom_node_point(aes(color = type), size = 5) +  # Color nodes by type
    geom_node_text(aes(label = name), repel = TRUE) +  # Add labels to nodes
    geom_edge_link(aes(edge_alpha = weight,edge_width = weight), color = "darkgrey") +  # Use edge weights for transparency
    theme_void() +  # Remove background and axes
    labs(title = "") +
    scale_color_manual(values = c("RS" = berryFunctions::addAlpha("blue",0.5),
                                  "US" = berryFunctions::addAlpha("orange",0.5),
                                  "DS" = berryFunctions::addAlpha("red",0.5),
                                  "SKIN" = berryFunctions::addAlpha("yellow",0.5),
                                  "Description_GO" = "lightgrey")) +  # Customize colors
    scale_edge_alpha(range = c(0.2, 1))+  # Adjust edge transparency based on weight
    scale_edge_width(range = c(0.3, 1.5))  # Adjust edge transparency based on weight
}

# Function to convert DataFrame to Cytoscape.js JSON
convert_to_cytoscape_json <- function(df) {

  library(jsonlite)

  # Create nodes and edges lists
  nodes <- list()
  edges <- list()

  # Add unique SysCT as a node
  unique_sysct <- unique(df$SysCT)
  for (sysct in unique_sysct) {
    nodes <- append(nodes, list(list(data = list(id = sysct, label = sysct))))
  }

  # Add unique Descriptions as nodes
  unique_descriptions <- unique(df$Description)
  for (desc in unique_descriptions) {
    nodes <- append(nodes, list(list(data = list(id = desc, label = desc))))
  }

  # Create edges based on SysCT and Description
  for (i in 1:nrow(df)) {
    edges <- append(edges, list(list(data = list(source = df$SysCT[i], target = df$Description[i]))))
  }

  # Create the final JSON structure
  cytoscape_json <- list(nodes = nodes, edges = edges)

  # Convert to JSON format
  json_output <- toJSON(cytoscape_json, pretty = TRUE, auto_unbox = TRUE)

  return(json_output)
}

### 202507 sys epi CNTN4 and SysCT correlation
subsysList=readRDS("./subsysList_AddGOscore.250220.rds")
inputV=c("GO:0016055", "GO:0001701", "GO:0019827", "GO:0048863", "GO:0030154", "GO:0048699", "GO:0048666", "GO:0030425", "GO:0007411")

inV=x01v
inV_name=x01v_name
skinV00=inV[inV %in% inputV]
skinV=inV_name[skinV00]
names(skinV)=names(skinV00)

for(x in names(skinV)){
#xxx=FeatureScatter(subsysList$SKIN, feature1 = x, feature2 = "CNTN4", group.by="refine_subtypes01")
#xxx=FeatureScatter(subsysList$DS, feature1 = x, feature2 = "CNTN4", group.by="refine_subtypes01")
#xxx=FeatureScatter(subsysList$US, feature1 = x, feature2 = "CNTN4", group.by="refine_subtypes01")
xxx=FeatureScatter(subsysList$RS, feature1 = x, feature2 = "CNTN4", group.by="refine_subtypes01")

xxx_data=xxx$data
xxx_data$Sys=gsub("_.+","",colnames(xxx_data)[1])
xxx_data$GO=skinV[x] %>% as.vector()
colnames(xxx_data)=c("Score","CNTN4","CT","Sys","GO")
allDF=rbind(allDF, xxx_data)
}

rownames(allDF)=NULL
ggplot()+
  geom_point(data=allDF[allDF$CT %in% c("CNTN4_NEL_Cells"),], aes(Score, CNTN4), alpha=0.4, size=1, color="darkred")+
  facet_wrap(~GO, ncol = 2)
