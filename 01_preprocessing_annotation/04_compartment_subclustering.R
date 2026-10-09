# ============================================================================
# Second-round sub-clustering of stromal / endothelial / lymphoid / myeloid compartments
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L4644-4957, L6892-7282)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
# NOTE: the original session sourced a personal package-loading helper here
#       (library(...) calls for Seurat, ggplot2, ComplexHeatmap, pheatmap, etc.).
#       Load the packages listed in the repository README instead.
setwd("./data/")
library(scPioneer)
library(reticulate)
reticulate::use_condaenv("SCS")

blank= theme(axis.title=element_blank(),
             axis.text=element_blank(),
             axis.ticks=element_blank())
myfont= theme(text=element_text(size=14,  family="Helvetica", colour ="black"))
NoTitle=theme(plot.title = element_blank())
ThemeBox=theme_bw()+theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

h5ad2seu=function(inh5, dataname="ABC", savedir="./"){

  sc <- import("scanpy")
  adata <- sc$read_h5ad(inh5)
  print(adata)

  srt_meta=as.data.frame(adata$obs)
  rownames(srt_meta)=paste0("cb_",rownames(srt_meta))  #srt_meta$cellid #
  head(srt_meta)

  Radata=adata$raw
  xmat <- t(Radata$X)
  if (!inherits(xmat, "dgCMatrix")) {
    xmat <- as.sparse(xmat[1:nrow(xmat), , drop = FALSE])
  }
  rownames(xmat) <- Radata$var_names$values %>% as.character()
  colnames(xmat) <- paste0("cb_", Radata$obs_names$values) #lapply(adata$obs_names$values, FUN = let2nums, l2n = FALSE) %>% unlist()
  srt <- CreateSeuratObject(counts = xmat, meta.data = srt_meta)
  srt@assays$RNA$data=srt@assays$RNA$counts

  for (x in c("X_umap","X_pca","X_pca_harmony")){
    obsm=adata$obsm[[x]]
    prefix_tmp=gsub("X_","",x)
    colnames(obsm) <- paste0(prefix_tmp, "_", seq_len(ncol(obsm)))
    rownames(obsm) <- paste0("cb_", adata$obs_names$values)
    srt[[prefix_tmp]] <- CreateDimReducObject(embeddings = obsm, assay = "RNA", key = prefix_tmp)
  }

  saveRDS(srt, file = paste0(savedir,dataname,"_seuratAnno.v0525.rds"))
  return(srt)

}

################################################################################

OutputDir="./AnaRes20240718/"
dir.create(OutputDir, recursive = T)

resDir=OutputDir
dataoi=c("strom", "endo", "lym", "myl")

df_DEGs=data.frame()
df_GOs=data.frame()
df_CellNums=data.frame()
for(i in 1:length(dataoi)){

  mydata=dataoi[i]

typecol="SubClass00" #paste0(shortCToi[i],"subtype00")
tmpseu=subdataList[[mydata]]
DefaultAssay(tmpseu)="RNA"
Idents(tmpseu)=typecol

deg000=DEG2FUN_subcts_list[[mydata]]$RawDEG[
  abs(DEG2FUN_subcts_list[[mydata]]$RawDEG$avg_log2FC)>0.25 &
    DEG2FUN_subcts_list[[mydata]]$RawDEG$p_val_adj<0.01,]
fun000=DEG2FUN_subcts_list[[mydata]]$FUNdf$GO

deg000$dataclass=mydata
df_DEGs=rbind(df_DEGs, deg000)
fun000$dataclass=mydata
df_GOs=rbind(df_GOs, fun000)

tmpseu$Location_short=droplevels(tmpseu$Location_short)
cnum000=table(tmpseu$SubClass00, tmpseu$Location_short) %>% as.data.frame()
cnum000$dataclass=mydata
df_CellNums=rbind(df_CellNums, cnum000)

}
subdataList$strom=subset(subdataList$strom, SubClass00 %!in% c("rmv"))
subdataList$strom$SubClass00=droplevels(subdataList$strom$SubClass00)
subdataList$lym=subset(subdataList$lym, SubClass00 %!in% c("CD8T_lq"))
subdataList$lym$SubClass00=droplevels(subdataList$lym$SubClass00)
subdataList$endo=subset(subdataList$endo, SubClass00 %!in% c("rmv"))
subdataList$endo$SubClass00=droplevels(subdataList$endo$SubClass00)
subdataList$strom=fastRefine(subdataList$strom, dims = 1:50, resol = 0.8, nnei = 100, spread = 1.2)
subdataList$endo=fastRefine(subdataList$endo, dims = 1:50, resol = 0.8, nnei = 100, spread = 1.2)

DimPlot(subdataList$endo, group.by = "SubClass00", reduction = "umap",label = T, repel = T)+NoLegend()+
  DimPlot(subdataList$endo, group.by = "Location_short", reduction = "umap",label = T)

Allgenes=rownames(subdataList$strom_fb@assays$RNA$data)
hoxgenes=Allgenes[grepl("^HOX", Allgenes)]
inobj=subset(subdataList$strom_fb, features=hoxgenes)
tmp_inobj=AverageExpression(
  inobj,
  assays = "RNA", group.by = "Location_short",return.seurat = T)
tmp_mtx=tmp_inobj@assays$RNA$data
tmp_mtx=tmp_mtx[rownames(tmp_mtx)[rowSums(tmp_mtx)>0.2],
                colnames(tmp_mtx)[colSums(tmp_mtx)>0.2]]
tmp_mtx=scale(tmp_mtx, center = TRUE, scale = TRUE) %>% as.matrix()
tmp_mtx=rescale(tmp_mtx, to = c(-2,2))

pheatmap(tmp_mtx, scale = TRUE, name="Z-score",
         cluster_rows = TRUE, cluster_cols = FALSE, )

################################################################################
inobj=subdataList$myl_mp
GO_DATA=clusterProfiler:::get_GO_data(org.Hs.eg.db, "BP", "SYMBOL")
SigName=c("ImmResp.","Chemok.")
SigGO=c("GO:0006955","GO:0070098")
for(x in 1:length(SigName)){
  inobj=AddModuleScore(inobj, features = list(GO_DATA$PATHID2EXTID[[SigGO[x]]]),
                      ctrl = 100, name = SigName[x])
}
M1_human <- c("NOS2", "CD86", "CD80", "CD38", "FCGR3A", "FCGR3B", "FCGR2B", "IL1B", "IL6", "IL12A", "IL12B", "IL23A", "TNF")
M2_human <- c("CD163", "MRC1", "IL10", "ARG1", "CHI3L1", "RETNLA", "TFRC", "TGFB1", "VEGFA", "SLC40A1", "NFE2L2", "HMOX1")
inobj=AddModuleScore(inobj, features = list(M1_human),
                     ctrl = 100, name = "M1.0")
inobj=AddModuleScore(inobj, features = list(M2_human),
                     ctrl = 100, name = "M2.0")

OutputDir="./AnaRes20241011/"
ClassActMP=c("IL1A", "IL1B", "CXCL8", "TNF", "IL23A","IL6","CCL3","CCL5","CXCL10","CCL23","MRC1")
AlterActMP=c("CCL22","CCL18","CCL17","IL1RN")
inobj=AddModuleScore(inobj, features = list(ClassActMP),
                     ctrl = 100, name = "ClassActMP.0")
inobj=AddModuleScore(inobj, features = list(AlterActMP),
                     ctrl = 100, name = "AlterActMP.0")

head(inobj@meta.data)
allpathname=paste(c("ClassActMP.0", "AlterActMP.0"),"1",sep="")
for(x in 1:length(allpathname)){
  xxx=VlnPlot(
              subset(inobj, MYL_SubClass04_bak %!in% c("rmv_BS","rmv_Intestine","rmv_skin","rmv_spleen","rmv_UDT")),
              features=c(allpathname[x]),
              group.by = "MYL_SubClass04_bak")+NoLegend()+xlab("")
  xxx[[1]]$layers[[3]]=xxx[[1]]$layers[[1]]
  xxx[[1]]$layers[[1]]=xxx[[1]]$layers[[2]]
  xxx[[1]]$layers[[2]]=xxx[[1]]$layers[[3]]
  xxx[[1]]$layers[[3]]=NULL
  xxx=xxx+geom_boxplot(color="white", width=0.3, outlier.alpha = 0)
  runFig(paste0(OutputDir,"Vlnplot.MPsub.",allpathname[x],".pathway_score"), xxx, 15, 3.75)
}

wgcnaRun=function(sc, prefix="tmp", fraction=0.1, submode=TRUE, cmplist=NULL, harmony=NULL, min_cells = 50,
                  ctname="seurat_clusters", ct_vec=NULL, gpname="Sample", gp_vec=NULL, baseDir="./", plot_mode=FALSE){

  library(devtools)
  library(GeneOverlap)
  library(UCell)
  library(ggrepel)
  library(igraph)
  library(ggforestplot)
  library(WGCNA)
  library(hdWGCNA)
  library(Seurat)
  library(cowplot)
  library(patchwork)

  dir.create(baseDir, recursive = T)

  set.seed(12345)
  theme_set(theme_cowplot())

  ### subtypes01 majortypes Sample

  GmodRes=list()

  harmonized=FALSE
  if(!is.null(harmony)){
    harmonized=TRUE
  }

  if(!plot_mode){

  if(!is.null(ct_vec)){
    sc=sc[,sc@meta.data[,ctname] %in% ct_vec]
    sc@meta.data[,ctname]=droplevels(sc@meta.data[,ctname])
  }
  if(!is.null(gp_vec)){
    sc=sc[,sc@meta.data[,gpname] %in% gp_vec]
    sc@meta.data[,gpname]=droplevels(sc@meta.data[,gpname])
  }

  print(table(sc@meta.data[,ctname],sc@meta.data[,gpname]))

  if(fraction>0){
    sc <- SetupForWGCNA(
      sc,
      gene_select = "fraction",
      fraction = fraction, # fraction of cells that a gene needs to be expressed in order to be included
      group.by=ctname,
      wgcna_name = paste0(prefix,"_hdWGCNA") # the name of the hdWGCNA experiment
    )
  }else{
    if(submode){
      sc=FindVariableFeatures(sc, method = "vst", nfeatures=2000)
    }
    sc <- SetupForWGCNA(
      sc,
      gene_select = "variable", # the gene selection approach
      wgcna_name = paste0(prefix,"_hdWGCNA") # the name of the hdWGCNA experiment
    )
  }
  reduc_name="pca"
  if(!is.null(harmony)){
    reduc_name=harmony
  }
  sc <- MetacellsByGroups(
    seurat_obj = sc,
    group.by = c(ctname, gpname), # specify the columns in seurat_obj@meta.data to group by
    reduction = reduc_name, # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    min_cells = min_cells,
    ident.group = ctname # set the Idents of the metacell seurat object
  )
  sc <- NormalizeMetacells(sc)
  #####################################################

  ct_vec=table(sc@misc[[paste0(prefix,"_hdWGCNA")]]$wgcna_metacell_obj[[ctname]][[1]]) %>% names()
  gp_vec=table(sc@misc[[paste0(prefix,"_hdWGCNA")]]$wgcna_metacell_obj[[gpname]][[1]]) %>% names()
  sc <- SetDatExpr(
    sc,
    group_name = ct_vec, # the name of the group of interest in the group.by column
    group.by=ctname, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data', # using normalized data
    multi.group.by = gpname,
    multi_group_name = gp_vec
  )

  sc <- TestSoftPowers(
    sc,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )

  sc <- ConstructNetwork(
    sc,
    overwrite_tom = TRUE,
    tom_outdir = baseDir,
    tom_name = prefix # name of the topoligical overlap matrix written to disk
  )

  # need to run ScaleData first or else harmony throws an error:

  # compute all MEs in the full single-cell dataset
  harmony_var=NULL
  if(!is.null(harmony)){
    harmony_var="Sample"
  }
  sc <- ModuleEigengenes(
    sc,
    group.by.vars=harmony_var # No harmony or #"Sample"
  )

  sc <- ModuleConnectivity(
    sc,
    group.by = ctname, harmonized = harmonized, group_name = ct_vec,
  )
  sc <- ResetModuleNames(
    sc,
    new_name = paste0(prefix,"-M")
  )

  #compute gene scoring for the top 25 hub genes by kME for each module
  sc <- ModuleExprScore(
    sc,
    n_genes = 25,
    method='UCell'
  )
  # get hMEs from seurat object
  t_MEs <- GetMEs(sc, harmonized=harmonized)
  # add hMEs to Seurat meta-data:
  sc@meta.data <- cbind(sc@meta.data, t_MEs)

  saveRDS(sc, file=paste0(baseDir,prefix,'.hdWGCNA_object.rds'))

  }
  ################################################################################

tmpTOM <- GetTOM(sc)
GmodRes[[paste0("TOM_matrix")]]=tmpTOM

MEs <- GetMEs(sc, harmonized=harmonized) #hMEs <- GetMEs(seurat_obj)
MEs_names=sort(colnames(MEs)); print(MEs_names)
GmodRes[[paste0("MEs_df")]]=MEs

modules <- GetModules(sc) %>% subset(module != 'grey') ###greys SHOULD be removed in ALL downstream analysis
# show the first 6 columns:
head(modules)
GmodRes[[paste0("modules_df")]]=modules

# get top25 hub genes
hub_df <- GetHubGenes(sc, n_hubs = 100)
head(hub_df)
GmodRes[[paste0("hub100genes_df")]]=hub_df

hub_df$gene=hub_df$gene_name
hub_df$avg_log2fc=100
hub_df$cluster=hub_df$module
ResFun=enrichFUNgo(ud_table = hub_df[hub_df$cluster %!in% c("grey"),],sp="hsa",intermode = FALSE, enriches = c("GO"))
GmodRes[[paste0("MEfun_GO")]]=ResFun$dflist$GO
GOdf=ResFun$dflist$GO
ddd=GOdf%>% group_by(celltype) %>% top_n(n=5, wt=-pvalue)
ddd$celltype=factor(ddd$celltype, levels = MEs_names[MEs_names %!in% c("grey")])
ddd$Description=factor(ddd$Description, levels = rev(unique(ddd$Description)))

GmodRes[[paste0("MEfun_GOplot")]]=ggplot(ddd)+
  geom_point(aes(celltype, Description, color=-log10(pvalue), size=Factor))+
  xlab("")+ylab("")+theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(paste0(baseDir,prefix,"_hdWGCNA_Dendrogram.pdf"), width = 8, height = 4)
PlotDendrogram(sc, main=paste0(prefix,' hdWGCNA Dendrogram'))
dev.off()

GmodRes[["KMEsPlot"]]=PlotKMEs(sc, ncol=4)

GmodRes[["ModuleRadarPlot"]]=ModuleRadarPlot(
  sc[,sc@meta.data[, ctname] %in% ct_vec],
  group.by = ctname,
  barcodes = sc[,sc@meta.data[, ctname] %in% ct_vec] %>% colnames(),
  axis.label.size=4,
  grid.label.size=4,
  combine=FALSE
)

## make a featureplot of hMEs for each module
#)
## stitch together with patchwork

# plot with Seurat's DotPlot function
t_modules <- GetModules(sc)
t_mods <- levels(t_modules$module); t_mods <- t_mods[t_mods != 'grey']
#xxx=DotPlot(sc, features=t_mods, group.by = ctname)+
GmodRes[["ModuleDotPlot_split"]]=splitDotPlot(sc[,sc@meta.data[, ctname] %in% ct_vec],
                                              group.by=ctname, groupV=NULL, split.by=gpname, splitV=gp_vec, featureV=t_mods)

  if(!is.null(cmplist)){

    # list of clusters to loop through
    clusters <- ct_vec
    cmplist=cmplist #list(c("B-A-3D","B"))
    # set up an empty dataframe for the DMEs
    DMEs <- data.frame()
    # all clusters as one
    for(i in 1:length(cmplist)){
      cmpitem=paste(cmplist[[i]], collapse = " vs. ", sep="")

      # identify barcodes for group1 and group2 in eadh cluster
      group1 <- sc[,sc@meta.data[, gpname] %in% cmplist[[i]][1] ] %>% colnames()
      group2 <- sc[,sc@meta.data[, gpname] %in% cmplist[[i]][2] ] %>% colnames()

      # run the DME test
      cur_DMEs <- FindDMEs(
        sc,
        barcodes1 = group1,
        barcodes2 = group2,
        test.use='wilcox',
        pseudocount.use=0.01, # we can also change the pseudocount with this param
        wgcna_name = sc@misc$active_wgcna #'INH'
      )

      # add the cluster info to the table
      cur_DMEs$cluster <- "All"
      cur_DMEs$cmp <- cmpitem

      # append the table
      DMEs <- rbind(DMEs, cur_DMEs)

    }
    # loop through the clusters
    for(cur_cluster in clusters){

      for(i in 1:length(cmplist)){
        cmpitem=paste(cmplist[[i]], collapse = " vs. ", sep="")

        # identify barcodes for group1 and group2 in eadh cluster
        group1 <- sc[,sc@meta.data[, ctname] %in% cur_cluster &
                       sc@meta.data[, gpname] %in% cmplist[[i]][1] ] %>% colnames()
        group2 <- sc[,sc@meta.data[, ctname] %in% cur_cluster &
                       sc@meta.data[, gpname] %in% cmplist[[i]][2] ] %>% colnames()

        if(length(group1)<5 | length(group2)<5) next

        # run the DME test
        cur_DMEs <- FindDMEs(
          sc,
          barcodes1 = group1,
          barcodes2 = group2,
          test.use='wilcox',
          pseudocount.use=0.01, # we can also change the pseudocount with this param
          wgcna_name = sc@misc$active_wgcna #'INH'
        )

        # add the cluster info to the table
        cur_DMEs$cluster <- cur_cluster
        cur_DMEs$cmp <- cmpitem

        # append the table
        DMEs <- rbind(DMEs, cur_DMEs)
      }
    }

    GmodRes[["InterPlot_split"]]=PlotDMEsLollipop(
      sc, #[,sc@meta.data[,ctname] %in% unique(DMEs$cluster)],
      DMEs,
      wgcna_name=sc@misc$active_wgcna,
      group.by = "cluster",
      comparison = c("All",ct_vec[ct_vec %in% unique(DMEs$cluster)]),
      pvalue = "p_val_adj"
    )

  }

  return(GmodRes)

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
