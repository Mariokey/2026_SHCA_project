# ============================================================================
# Cell-type specific location-similarity dendrograms (top HVGs, Euclidean distance, Ward clustering)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Dendrogram comparison"
# Source     : extracted from Analysis240525.r (L2237-2253, L12495-12557)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### cell_hierarchy tree

library(Census)
ggsave(plot_cell_hierarchy(Census:::census_ts_model$hierarchy_mat,
                           axis.limits = c(15,-15),
                           edge.thickness = 0.5,
                           label.line.padding = 0.1,
                           text.size = 1.5,
                           label.size = 0.75
),
filename = 'cell_hierarchy.pdf', width = 3, height = 8)

res = census_main(subseuList$Lymphocytes, organ = NULL, controu_data=TRUE)

### 20250812

x=43
tmpdata=subset(Mobj_keep, subtypes015 %in% levels(Mobj_keep$subtypes015)[x])
#}
tmpdata$Location_short=droplevels(tmpdata$Location_short)
tmpdata=FindVariableFeatures(tmpdata)
LS_Gene <- VariableFeatures(object = tmpdata)
tdata <- subset(tmpdata,features = LS_Gene)
Idents(tdata)="Location_short"
tdata <- AverageExpression(tdata,
                           return.seurat = TRUE)
Mat_Tmp <- as.matrix(tdata@assays$RNA$data)

ct_color=hue_pal()(dim(Mat_Tmp)[2])
names(ct_color)=colnames(Mat_Tmp)

hc= hclust(dist(t(Mat_Tmp)))
p <- ggtree(hc, linetype='solid')
d <- data.frame(label = colnames(Mat_Tmp),
                cyl =colnames(Mat_Tmp))

p = p %<+% d +
  layout_dendrogram() +
  geom_tippoint(aes(fill=factor(cyl), x=x+.5),
                size=8, shape=21, color='white') +
  geom_tiplab(angle=45, hjust=1, offset=-3, show.legend=FALSE) +
  scale_fill_manual(values = ct_color)+
  theme_tree(plot.margin=margin(50,50,80,50))+
  theme(legend.position="none")

reftree_df=read.xlsx("./organ_combinations_with_relationship_updated.xlsx")
reftree_df$organ1=gsub("Gallbladder","GBD",reftree_df$organ1)
reftree_df$organ2=gsub("Gallbladder","GBD",reftree_df$organ2)
reftree_df$organ1=gsub("\\s+","",reftree_df$organ1)
reftree_df$organ2=gsub("\\s+","",reftree_df$organ2)
reftree=tidyr::pivot_wider(reftree_df, id_cols = "organ1", names_from  = "organ2", values_from = "relationship")
reftree_mat=reftree[,2:dim(reftree)[2]] %>% as.matrix()
rownames(reftree_mat)=reftree[[1]] %>% as.vector()
xxx=pheatmap(reftree_mat, name = "ref_dis", clustering_method = "wald.D")
