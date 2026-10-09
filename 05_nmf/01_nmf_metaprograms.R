# ============================================================================
# NMF factorisation, rank survey, meta-program calling and CNTN4+ program analysis
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "NMF analysis"
# Source     : extracted from Analysis240525.r (L8950-9202, L10544-11046)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### NMF with filter data
xxx <- NormalizeData(xxx, normalization.method = "LogNormalize", scale.factor = 10000)
xxx <- FindVariableFeatures(xxx, selection.method = "vst", nfeatures = 3000)
xxx <- ScaleData(object = xxx, vars.to.regress = "pct_counts_mt")
hvg <- VariableFeatures(xxx, selection.method = "vst")[1:3000]

Idents(xxx) <- xxx$Location_short
keep_set <- as.data.frame(table(xxx$Location_short)) %>% filter(Freq > 50)

library(NMF)  ###Try loading NMF package before running
nmf_out <- list()
for( i in 1:length(keep_set$Var1)){
  patient = keep_set$Var1[i]
  print(patient)
  Subset <- subset(xxx, idents=patient)
  print(dim(Subset@meta.data))
  Subset <- NormalizeData(Subset, normalization.method = "LogNormalize", scale.factor = 10000)
  Subset <- ScaleData(object = Subset, features = hvg, vars.to.regress = "pct_counts_mt")
  tdata <- as.matrix(Subset@assays$RNA$scale.data)
  tdata[tdata<0] <- 0
  tdata <- tdata[rowSums(tdata) > 0, ]
  set.seed(123)
  NMFs_per_sample = NMF::nmf(x = tdata, rank = 4:9, method="snmf/r", nrun = 10)
  save(NMFs_per_sample, file = paste0("./AnaRes20241220/",patient,"_SP_NMF.RData"))
  Genes_nmf_w_basis <- cbind(NMFs_per_sample$fit$`4`@fit@W,
                             NMFs_per_sample$fit$`5`@fit@W,
                             NMFs_per_sample$fit$`6`@fit@W,
                             NMFs_per_sample$fit$`7`@fit@W,
                             NMFs_per_sample$fit$`8`@fit@W,
                             NMFs_per_sample$fit$`9`@fit@W)

  colnames(Genes_nmf_w_basis) <- c(paste0(patient,"_rank4_", seq(1:4)),
                                   paste0(patient,"_rank5_", seq(1:5)),
                                   paste0(patient,"_rank6_", seq(1:6)),
                                   paste0(patient,"_rank7_", seq(1:7)),
                                   paste0(patient,"_rank8_", seq(1:8)),
                                   paste0(patient,"_rank9_", seq(1:9)))

  nmf_out[[i]] <- Genes_nmf_w_basis

}
save(nmf_out, file = "./AnaRes20241220/Genes_nmf_All.rds")
### postprocess
source("./NMF/20240000_custom_magma.R")
source("./NMF/20240000_robust_nmf_programs.R")

## Parameters
intra_min_parameter <- 35
intra_max_parameter <- 10
inter_min_parameter <- 10
names(nmf_out) <- keep_set$Var1   ####
nmf_programs  <- lapply(nmf_out, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))

nmf_programs <- lapply(nmf_programs,toupper) ## convert all genes to uppercase
nmf_filter_ccle <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

# calculate similarity between programs
nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y))))

# hierarchical clustering of the similarity matrix
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average")
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------

Genes_nmf_w_basis <- nmf_out

### Parameters for clustering
Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 10    # the minimal intersection cutoff for adding a new NMF to the forming cluster
Min_group_size        <- 5     # the minimal group size to consider for defining the first NMF program in a cluster

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {

  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))

  ### intersection between all remaining NMFs and Genes in MP
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF

  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {

    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])

    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene

    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste(strsplit(i , "_rank")[[1]][1])
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i]
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]

      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))

    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50]
    }

    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])
    Genes_MP        <- Genes_MP_temp[1:50]

    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF

    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP

  }

  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1

  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster

  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters

  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

####  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){

  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))

}

inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new])
library(scales)
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
ggsave("./AnaRes20241220/byLocation_Epi_cells-10MetaProgram.pdf", width = 8, height = 7)

MP_list <-  do.call(cbind, MP_list)

MP_genes_DF=MP_list %>% as.data.frame() %>% tidyr::pivot_longer(cols = MP_1:MP_12)
colnames(MP_genes_DF)=c("cluster","gene")
MP_genes_DF$comparison=MP_genes_DF$cluster
MP_genes_DF$UPfor=MP_genes_DF$cluster
MP_genes_DF$avg_log2fc=10
MP_funs_list=enrichFUNgo(MP_genes_DF, intermode = FALSE, enriches = c("GO","KEGG"))

GOfun=MP_funs_list$dflist$KEGG %>% group_by(celltype) %>% top_n(n=10, wt=-pvalue) %>% as.data.frame()
GOfun$Factor=sapply(GOfun$GeneRatio, function(x) eval(parse(text=x)))
ggplot(GOfun)+geom_point(aes(Factor, Descript, size=Count, color=pvalue))+
  scale_color_gradient(low="red", high="blue", trans="log10")+
  facet_wrap(~celltype, ncol = 4, scales = "free")+
  xlab("Rich Factor")+ylab("KEGG Term")

for(x in 1:length(unique(MP_genes_DF$cluster))){
  MPid=unique(MP_genes_DF$cluster)[x]
  genes_vec=MP_genes_DF[MP_genes_DF$cluster %in% MPid,]$gene
  seu_epi_curat=AddModuleScore(seu_epi_curat, assay = "RNA", name = paste0(MPid,"."),
                               features = list(genes_vec), ctrl = 100)
}

MP12score_DF=seu_epi_curat@meta.data[,c("forCombine01","Location_short",paste0(unique(MP_genes_DF$cluster),".1"))]
MP12score_DF_longer=MP12score_DF %>% tidyr::pivot_longer(cols=MP_1.1:MP_12.1) %>% as.data.frame()
MP12score_DF_longer_mean=MP12score_DF_longer %>% group_by(forCombine01,Location_short,name) %>% summarise_at("value",mean)

MP12score_DF_longer_mean=MP12score_DF_longer_mean[order(MP12score_DF_longer_mean$forCombine01,MP12score_DF_longer_mean$Location_short,MP12score_DF_longer_mean$value, decreasing = c(FALSE, FALSE, TRUE)),]
ggplot(MP12score_DF_longer_mean)+geom_point(aes(name, value))+
  facet_grid(Location_short~forCombine01, scales = "free")+
  theme(axis.title = element_text(size=14, colour = "black"),
        axis.text = element_text(size=14, colour = "black"),
        strip.text = element_text(size=14, colour = "black"),
        strip.background = element_blank())

MP12score_DF_longer_mean$Location_short=factor(
  MP12score_DF_longer_mean$Location_short, levels = levels(seu_epi_curat$Location_short)
)
MP12score_DF_longer_mean$forCombine01=factor(
  MP12score_DF_longer_mean$forCombine01, levels = levels(seu_epi_curat$forCombine01)
)

plot_tmplist=list()
for(x in unique(MP12score_DF_longer_mean$name)){
  tmp_mean=MP12score_DF_longer_mean[MP12score_DF_longer_mean$name %in% c(x),]
  print(table(tmp_mean$name))
  pp0=ggplot(data = tmp_mean[tmp_mean$value>0,], aes(x=Location_short, y=forCombine01, fill=value, color=value)) +
    geom_tile() +
    scale_color_gradient2(limits=c(min(tmp_mean$value),max(tmp_mean$value)), low=custom_magma[333:223],  mid =custom_magma[222:112], high = custom_magma[111:1], midpoint = 0, oob=squish, name=paste0(gsub("\\.1","",x),"\nScore")) +
    scale_fill_gradient2(limits=c(min(tmp_mean$value),max(tmp_mean$value)), low=custom_magma[333:223],  mid =custom_magma[222:112], high = custom_magma[111:1], midpoint = 0, oob=squish, name=paste0(gsub("\\.1","",x),"\nScore"))  +
    theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, color="black"), axis.text.y=element_text(color="black"))+
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
  print(pp0)
  ggsave(paste0("./AnaRes20241220/",x,"_LocByCT.pdf"), width = 18, height = 8)
  plot_tmplist[[x]]=pp0
}

## CNTN4+
CNTN4_data=fastRefine(CNTN4_data, dims = 1:50, resol = 0.8, nnei = 100, spread = 1.2)
DimPlot(CNTN4_data, group.by = "seurat_clusters", reduction = "umap",label = T)
scPioneer::FeaturePlot2(CNTN4_data, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
scPioneer::DimPlot_idx(CNTN4_data,  group.by = "seurat_clusters", reduction = "umap")

CNTN4_data$System=factor(CNTN4_data$System, levels = levels(Mobj_keep$System))
CNTN4_data$Location_short=factor(CNTN4_data$Location_short, levels = levels(Mobj_keep$Location_short))
inobj=CNTN4_data
inobj$UMAP_1=inobj@reductions$umap@cell.embeddings[,1]
inobj$UMAP_2=inobj@reductions$umap@cell.embeddings[,2]
gpColor=hue_pal()(length(levels(inobj$System)))
names(gpColor)=levels(inobj$System)
plotlist=list()
for(x in names(gpColor)){
  ppp=ggplot()+stat_density_2d(data=inobj@meta.data, #%>% sample_n(size=10000),
                               aes(UMAP_1, UMAP_2, color=as.factor("..level..")), breaks=c(0.0005,100), h = c(0.8, 0.8), n=100)+
    scale_color_manual(values = c("grey",rep("NA",100),"NA"))+
    stat_density_2d(data=inobj@meta.data[inobj@meta.data$System %in% x, ],
                    aes_string(x = "UMAP_1", y = "UMAP_2",  alpha = "..level..", fill = "System"),
                    linewidth = 0, geom = "polygon",  n = 200, h = c(1, 1)) + theme_cowplot() +
    scale_fill_manual(values = gpColor)+
    scale_alpha_continuous(limits=c(0,0.3)) + NoLegend()
  plotlist[[x]]=ppp
}
plot_grid(plotlist = plotlist, ncol = 3, labels = names(plotlist))

inobj$seurat_clusters=as.numeric(as.character(inobj$seurat_clusters))
inobj$seurat_clusters=inobj$seurat_clusters+1
inobj$seurat_clusters=factor(inobj$seurat_clusters, levels = 1:10)
inobj$Location_short=droplevels(inobj$Location_short)
roveRun(inobj, Group="Location_short", CT="seurat_clusters", outdir=paste0(outdir,"/LocationCNTN4.clu10"), flip = FALSE, mywidth=10, nCT=10)

### NMF with filter data
newDir="./AnaRes20250215/"
dir.create(newDir)
xxx=inobj
xxx <- NormalizeData(xxx, normalization.method = "LogNormalize", scale.factor = 10000)
xxx <- FindVariableFeatures(xxx, selection.method = "vst", nfeatures = 3000)
xxx <- ScaleData(object = xxx, vars.to.regress = "pct_counts_mt")
hvg <- VariableFeatures(xxx, selection.method = "vst")[1:3000]

Idents(xxx) <- xxx$Location_short #xxx$Organ
keep_set <- as.data.frame(table(xxx$Location_short)) %>% filter(Freq > 10)

library(NMF)  ###Try loading NMF package before running
nmf_out <- list()
for( i in 1:length(keep_set$Var1)){
patient = keep_set$Var1[i]
print(patient)
Subset <- subset(xxx, idents=patient)
print(dim(Subset@meta.data))
Subset <- NormalizeData(Subset, normalization.method = "LogNormalize", scale.factor = 10000)
Subset <- ScaleData(object = Subset, features = hvg, vars.to.regress = "pct_counts_mt")
tdata <- as.matrix(Subset@assays$RNA$scale.data)
tdata[tdata<0] <- 0
tdata <- tdata[rowSums(tdata) > 0, ]
set.seed(123)
NMFs_per_sample = NMF::nmf(x = tdata, rank = 4:9, method="snmf/r", nrun = 10)
save(NMFs_per_sample, file = paste0(newDir,patient,"_SP_NMF.RData"))

Genes_nmf_w_basis=data.frame()
for(t_i in names(NMFs_per_sample$fit)){
  #t_i=4
  t_i=as.character(t_i)

  if(!is.na(NMFs_per_sample$fit[[t_i]])){

    t_g2nmf=NMFs_per_sample$fit[[t_i]]@fit@W
    colnames(t_g2nmf)=paste0(patient,"_rank",t_i,"_", seq(1:t_i))
    if(dim(Genes_nmf_w_basis)[2]==0){
      Genes_nmf_w_basis=t_g2nmf
    }else{
      Genes_nmf_w_basis=cbind(Genes_nmf_w_basis,t_g2nmf)
    }

  }

}

nmf_out[[i]] <- Genes_nmf_w_basis

}
save(nmf_out, file = paste0(newDir,"/Genes_nmf_All.rds"))
### postprocess
source("./NMF/20240000_custom_magma.R")
source("./NMF/20240000_robust_nmf_programs.R")

## Parameters
intra_min_parameter <- 35
intra_max_parameter <- 10
inter_min_parameter <- 10
names(nmf_out) <- keep_set$Var1   ####
nmf_programs  <- lapply(nmf_out, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))

nmf_programs <- lapply(nmf_programs,toupper) ## convert all genes to uppercase
nmf_filter_ccle <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

# calculate similarity between programs
nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y))))

# hierarchical clustering of the similarity matrix
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average")
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------

Genes_nmf_w_basis <- nmf_out

### Parameters for clustering
Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 10    # the minimal intersection cutoff for adding a new NMF to the forming cluster
Min_group_size        <- 5     # the minimal group size to consider for defining the first NMF program in a cluster

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {

  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))

  ### intersection between all remaining NMFs and Genes in MP
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF

  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {

    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])

    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene

    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste(strsplit(i , "_rank")[[1]][1])
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i]
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]

      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))

    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50]
    }

    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])
    Genes_MP        <- Genes_MP_temp[1:50]

    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF

    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP

  }

  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1

  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster

  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters

  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

####  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){

  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))

}

inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new])
library(scales)
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
ggsave(paste0(newDir,"/byOrganCNTN4_cells-10MetaProgram.pdf"), width = 8, height = 7)

MP_list <-  do.call(cbind, MP_list)

MP_genes_DF=MP_list %>% as.data.frame() %>% tidyr::pivot_longer(cols = MP_1:MP_5)
colnames(MP_genes_DF)=c("cluster","gene")
MP_genes_DF$comparison=MP_genes_DF$cluster
MP_genes_DF$UPfor=MP_genes_DF$cluster
MP_genes_DF$avg_log2fc=10
MP_funs_list=enrichFUNgo(MP_genes_DF, intermode = FALSE, enriches = c("GO","KEGG"))

MP_funs_list=CNTN4_NMF_results_list$MP_funs_list
GOfun=MP_funs_list$dflist$GO %>% group_by(celltype) %>% top_n(n=10, wt=-pvalue) %>% as.data.frame()
GOfun$Factor=sapply(GOfun$GeneRatio, function(x) eval(parse(text=x)))
ggplot(GOfun)+geom_point(aes(Factor, Descript, size=Count, color=pvalue))+
scale_color_gradient(low="red", high="blue", trans="log10")+
facet_wrap(~celltype, ncol = 5, scales = "free")+
xlab("Rich Factor")+ylab("GO Term")

CNTN4_NMF_results_list=readRDS("./AnaRes20250215/CNTN4_NMF_results_list.rds")
library(scales)
ggplot(data = CNTN4_NMF_results_list$nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
ggsave(paste0(newDir,"/byOrganCNTN4_cells-10MetaProgram_label.pdf"), width = 12, height = 10)

tmplist=lapply(1:5, function(x){CNTN4_NMF_results_list$MP_list[,x]})
names(tmplist)=colnames(CNTN4_NMF_results_list$MP_list)
CNTN4_data=AddModuleScore(CNTN4_data, features = tmplist,
                               ctrl = 100)
CNTN4_data$Location_short=factor(CNTN4_data$Location_short, levels = levels(Mobj_keep$Location_short))
VlnPlot(CNTN4_data, features = c(paste0("Cluster",1:5)), group.by = "Location_short", stack = T)+geom_boxplot()

inobj$Organ=inobj@meta.data[,"Location_short"]
for(i in 1:length(refine_list)){
  for(j in 1:length(refine_list[[i]])){
    inobj=modifyAno(inobj, anoCol = "Organ", oldID = paste0("^",refine_list[[i]][j],"$"), newID = names(refine_list[i]))
  }
}
inobj$Organ=factor(inobj$Organ, levels = rev(names(refine_list)))
inobj$Organ=droplevels(inobj$Organ)

dir.create("./hdWGCNA2/")
inobj=subset(inobj, Organ %!in% c("Spleen","Testis","Lung","Liver"))    ############## CNTN4_data_refine.rds without four organs (fewer than 30 cells)
inobj=subset(inobj, seurat_clusters %!in% c("10"))                      #### 44 cells, only in "VD" location
inobj=subset(inobj, seurat_clusters %!in% c("9"))                       #### 140 cells, almost only in "Renal" location
inobj=subset(inobj, seurat_clusters %!in% c("6"))                       #### 450 cells, almost only in "Ovary" location
inobj$Organ=droplevels(inobj$Organ)
inobj$seurat_clusters=droplevels(inobj$seurat_clusters)
wgcna_CNTN4=wgcnaRun(inobj, min_cells = 25,   ###knn k=25 default
                  baseDir="./hdWGCNA2/", harmony="pca_harmony",
                  prefix = "MP", fraction=0.01,
                  submode=TRUE, cmplist=NULL,
                  ctname="seurat_clusters",
                  ct_vec = NULL,
                  gpname="Organ")    #### wgcna run WRONG, cell number too small !!! give up...?

### note: use all CNTN4_data cells here !!!
CNTN4_results <- CytoTRACE::CytoTRACE(as.matrix(CNTN4_data@assays$RNA$counts),
                           ncores = 4, subsamplesize = 4000)
ptV=CNTN4_data$seurat_clusters %>% as.character()
names(ptV)=names(CNTN4_data$seurat_clusters)
### NOTE DO NOT load scPioneer!!! ->boxplot will be replaced by this package!!!
CytoTRACE::plotCytoTRACE(CNTN4_results, phenotype = ptV,  outputDir = newDir,
              emb = CNTN4_data@reductions$umap@cell.embeddings %>% as.data.frame()) #, gene = "Kit")

mn2.CNTN4_data=runMonocle2(CNTN4_data, "seurat_clusters", mygroup = "Location_short", mode = c("disp"),
                               rmvIdent=c("10"), savefix="cds_CNTN4_250215.")
mn2.CNTN4_data=readRDS("./CNTN4_data_mn2.rds")

### 20250304
CNTN4_data=readRDS("./CNTN4_data_refine.rds")
CNTN4_NMF_results_list=readRDS("./AnaRes20250215/CNTN4_NMF_results_list.rds")

tmplist=lapply(1:5, function(x){CNTN4_NMF_results_list$MP_list[,x]})
names(tmplist)=colnames(CNTN4_NMF_results_list$MP_list)
CNTN4_data=AddModuleScore(CNTN4_data, features = tmplist,
                          ctrl = 100)
CNTN4_data$Location_short=factor(CNTN4_data$Location_short, levels = levels(Mobj_keep$Location_short))
VlnPlot(CNTN4_data, features = c(paste0("Cluster",1:5)), group.by = "Location_short", stack = T)+geom_boxplot()
VlnPlot(CNTN4_data, features = c(paste0("Cluster",1:5)), group.by = "seurat_clusters", stack = T)+geom_boxplot()

mn2.CNTN4_data=readRDS("./CNTN4_data_mn2.rds")

### 202508
CNTN4_data_ori=readRDS("./CNTN4_data.rds")
CNTN4_data_ori=fastRefine(CNTN4_data_ori, dims = 1:50, resol = 0.8, nnei = 100, spread = 1.2)
CNTN4_data_ori$System=factor(CNTN4_data_ori$System, levels = levels(Mobj_keep$System))
CNTN4_data_ori$Location_short=factor(CNTN4_data_ori$Location_short, levels = levels(Mobj_keep$Location_short))

CNTN4_data_ori$seurat_clusters=as.numeric(as.character(CNTN4_data_ori$seurat_clusters))
CNTN4_data_ori$seurat_clusters=CNTN4_data_ori$seurat_clusters+1
CNTN4_data_ori$seurat_clusters=factor(CNTN4_data_ori$seurat_clusters, levels = 1:10)
CNTN4_data_ori$Location_short=droplevels(CNTN4_data_ori$Location_short)
CNTN4_data_ori=subset(CNTN4_data_ori,seurat_clusters %!in% c("10") )
CNTN4_data_ori$seurat_clusters=droplevels(CNTN4_data_ori$seurat_clusters)

CNTN4_data_ori$types=CNTN4_data_ori$seurat_clusters
CNTN4_data_ori$cellnames=CNTN4_data_ori$cellid
CNTN4_data_ori$cellid=rownames(CNTN4_data_ori@meta.data)
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(1,2))), newID = "NEL01_MECOM")
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(3,4))), newID = "NEL02_KRT14")
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(5,6))), newID = "NEL03_ATRNL1")
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(7))), newID = "NEL04_CNTNAP2")
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(8))), newID = "NEL05_S100A2")
CNTN4_data_ori=modifyAno(CNTN4_data_ori, anoCol = "types", cellid = colnames(subset(CNTN4_data_ori, seurat_clusters %in% c(9))), newID = "NEL06_FXYD2")
CNTN4_data_ori$types=factor(CNTN4_data_ori$types, levels = sort(unique(CNTN4_data_ori$types)))

roveRun(CNTN4_data_ori, Group="Location_short", CT="types", outdir=paste0(CCIdir_img,"/LocationCNTN4.clu6"), flip = FALSE, mywidth=10, nCT=6)
Idents(CNTN4_data_ori)="types"
CNTN4_data_ori_degs=FindAllMarkers(CNTN4_data_ori, group.by = "types", only.pos = T, max.cells.per.ident = 2000, test.use = "t")
write.xlsx(CNTN4_data_ori_degs, file=paste0(CCIdir_img, "/CNTN4_merge6_degs.xlsx"))
CNTN4_data_ori_degs=read.xlsx("./CCI_sub0310/CNTN4_merge6_degs.xlsx")
DEGcl_top=as.data.frame(CNTN4_data_ori_degs %>% group_by(cluster) %>% top_n(n=5, wt=(-p_val_adj)))
xxx=scPioneer::DotPlot4(CNTN4_data_ori, features = unique(DEGcl_top$gene), group.by = "types")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
runFig(paste0(FinalFig,"/", "CNTN4marker.DotPlot_DEG"), xxx, 12, 4.5)

saveRDS(CNTN4_data_ori, file = "./CNTN4_data_ori.rds")

### CNTN4 and neuroendo。。。20250720
cntn4_6ct=sample(scPioneer::scPalette2(length(levels(CNTN4_data_ori@meta.data[, "types"]))))
names(cntn4_6ct)=levels(CNTN4_data_ori@meta.data[, "types"])
p=scPioneer::DimPlot_idx(CNTN4_data_ori, group.by = "types", reduction = "umap", cols=as.vector(cntn4_6ct))+NoAxes()
runFig("./CNTN4_output/CNTN4CT6.dimplot", p, 5.25,3.25)

p=scPioneer::FeaturePlot2(CNTN4_data_ori, features = c("EPCAM","CDH1","KRT18","CNTN4","PDZRN3","YAP1","GPC6", "FBXL7"), ncol = 4, order = T)
runFig("./CNTN4_output/CNTN4CT6.featureplot", p, 18,5.5)

CNTN4_data_ori$System=factor(CNTN4_data_ori$System, levels = levels(Mobj_keep$System))
CNTN4_data_ori$Location_short=factor(CNTN4_data_ori$Location_short, levels = levels(Mobj_keep$Location_short))
inobj=CNTN4_data_ori
inobj$UMAP_1=inobj@reductions$umap@cell.embeddings[,1]
inobj$UMAP_2=inobj@reductions$umap@cell.embeddings[,2]
gpColor=hue_pal()(length(levels(inobj$System)))
names(gpColor)=levels(inobj$System)
plotlist=list()
for(x in names(gpColor)){
  ppp=ggplot()+stat_density_2d(data=inobj@meta.data, #%>% sample_n(size=10000),
                               aes(UMAP_1, UMAP_2, color=as.factor("..level..")), breaks=c(0.0005,100), h = c(0.8, 0.8), n=100)+
    scale_color_manual(values = c("grey",rep("NA",100),"NA"))+
    stat_density_2d(data=inobj@meta.data[inobj@meta.data$System %in% x, ],
                    aes_string(x = "UMAP_1", y = "UMAP_2",  alpha = "..level..", fill = "System"),
                    linewidth = 0, geom = "polygon",  n = 200, h = c(1, 1)) + theme_cowplot() +
    scale_fill_manual(values = gpColor)+
    scale_alpha_continuous(limits=c(0,0.3)) + NoLegend()
  plotlist[[x]]=ppp
}
p=plot_grid(plotlist = plotlist, ncol = 3, labels = names(plotlist))
runFig("./CNTN4_output/CNTN4CT6.sysDistri_plot", p, 10,5.5)

inobj$Location_short=droplevels(inobj$Location_short)
roveRun(inobj, Group="Location_short", CT="types", outdir="./CNTN4_output/CNTN4CT6.rove", flip = FALSE, mywidth=10, nCT=7)

### note: use all CNTN4_data cells here !!!
CNTN4_results <- CytoTRACE::CytoTRACE(as.matrix(CNTN4_data_ori@assays$RNA$counts),
                                      ncores = 4, subsamplesize = 2000)
ptV=CNTN4_data_ori$types %>% as.character()
names(ptV)=names(CNTN4_data_ori$types)
### NOTE DO NOT load scPioneer!!! ->boxplot will be replaced by this package!!!
CytoTRACE::plotCytoTRACE(CNTN4_results, phenotype = ptV,  #colors = cntn4_6ct,
                         outputDir = "./CNTN4_output/",
                         emb = CNTN4_data_ori@reductions$umap@cell.embeddings %>% as.data.frame()) #, gene = "Kit")

sub_CNTN4_data_ori=subset(CNTN4_data_ori, types %in% c("NEL01_MECOM","NEL04_CNTNAP2"))
sub_CNTN4_data_ori=subset(sub_CNTN4_data_ori,
                          Location_short %in% levels(CNTN4_data_ori$Location_short)[1:22])
sub_CNTN4_data_ori$Location_short=droplevels(sub_CNTN4_data_ori$Location_short)
sub_CNTN4_data_ori=RunUMAP(sub_CNTN4_data_ori,dims = 1:50, reduction="pca_harmony", reduction.name = "umap",
                           n.neighbors = 100, min.dist = 0.15, spread = 1.2, verbose = FALSE)
DimPlot(sub_CNTN4_data_ori, group.by = "types", reduction = "umap",label = T)
mn2.CNTN4_data_subDS=runMonocle2(sub_CNTN4_data_ori, "types", mygroup = "Location_short", mode = c("disp"),
                           rmvIdent=NULL, savefix="./CNTN4_output/cds_CNTN4_250830.")
mn2.CNTN4_data_subDS$disp$pseudotime$data$Pseudotime_bak=mn2.CNTN4_data_subDS$disp$pseudotime$data$Pseudotime
mn2.CNTN4_data_subDS$disp$pseudotime$data$Pseudotime=max(mn2.CNTN4_data_subDS$disp$pseudotime$data$Pseudotime)-mn2.CNTN4_data_subDS$disp$pseudotime$data$Pseudotime

runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS2types_state"),
       mn2.CNTN4_data_subDS$disp$state, 3.5, 3.25)
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS2types_time"),
       mn2.CNTN4_data_subDS$disp$pseudotime, 3.5, 3.25)
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS2types_celltypes"),
       mn2.CNTN4_data_subDS$disp$celltype, 3.5, 3.25)
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS2types_Heatmap_celltypesXstate"),
       mn2.CNTN4_data_subDS$disp$heat_comb, 1.5, 3.25)

sub_CNTN4_data_ori$Cellname=rownames(sub_CNTN4_data_ori@meta.data)
sub_CNTN4_data_ori_sls_res0901=slingRun(sub_CNTN4_data_ori, samplesize=1, CTname="types", cellid="Cellname",
                           start.type=c("NEL01_MECOM"),
                           mycolor89=cntn4_6ct, savesce=paste0("./CNTN4_output/CNTN4_DS2cts_slingshot.0901.rds"))
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS2types_slingshot"),
       sub_CNTN4_data_ori_sls_res0901$UMAP_allCurve+theme(panel.grid  = element_blank()), 4, 2.5)
sce=readRDS("./CNTN4_output/CNTN4_DS2cts_slingshot.0901.rds")
sub_CNTN4_data_ori_tradeResList=tradeS(sce, samplesize=1, SvE=FALSE,  Asso=TRUE, outdirprefix=paste0("./CNTN4_output/tradeSeq0901/CNTN4_DS2cts"))

p=gseaplot(sub_CNTN4_data_ori_tradeResList$L1_GSEres$gse_obj,
         by = "all", title = name, geneSetID = name_go)
runFig(paste0(OutputDir,"/", paste0("Fig3.CNTN4_subDS_tradeseq_",gsub(":","",name_go))), p, 4.75, 5)

tmpheat=sub_CNTN4_data_ori_tradeResList$L1_smoothExp
tflist=read.table("./ref/human.TF.list")
target_genes=intersect(rownames(tmpheat),tflist$V1)
target_genes=target_genes[target_genes %!in% c("ZNF131","POU6F2","NR5A2")]
tmp_tf=tmpheat[target_genes,]
p=pheatmap(tmp_tf, cluster_rows = F, cluster_cols = F, name="TF_exp.",
         show_rownames = T, show_colnames = F)
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS_TF"), draw(p), 4.5, 6.25)

tmpheat=sub_CNTN4_data_ori_tradeResList$L1_smoothExp_heatplot
runFig(paste0(OutputDir,"/", "Fig3.CNTN4_subDS_allGenes"), draw(tmpheat), 4.5, 6.25)
