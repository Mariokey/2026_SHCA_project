# ============================================================================
# basicDraw00: composite overview figure per cell type (UMAP, Ro/e, DEG heatmap, enrichment)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Statistics and reproducibility"
# Source     : extracted from Analysis240525.r (L776-1235)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
###****** *
###**** *
###** *
basicDraw00=function(inobj, figfix="", iCT="Major6CT", rmvCT=NULL,
                   degDF=NULL, funDF=NULL, idegDF=NULL, ifunDF=NULL,
                   outdir="./", pcgeneDF=NULL, CTorder=NULL, GPorder=NULL, ...){

library(entropy)

######################  ######################   ######################
Gorder_file=read.xlsx("./../bin/meta_table.uniqorder.0505.xlsx")

CTorder=NULL
GPorder=NULL
pcgeneDF=read.table("./ref/hum_pcgene.list")

DefaultAssay(inobj)="RNA"
Idents(inobj)=iCT
if(!is.null(rmvCT)){
  inobj=inobj[,(inobj@meta.data[,iCT] %!in% rmvCT)]
}
inobj@meta.data[,iCT]=droplevels(inobj@meta.data[,iCT])
nCT=length(levels(inobj@meta.data[,iCT]))
if(is.null(CTorder)){
  CTorder=levels(inobj@meta.data[,iCT])
}
inobj$Location_short=factor(inobj$Location_short, levels = unique(Gorder_file$Location_short))
inobj=subset(inobj, Location_short %in% names(table(inobj$Location_short)[table(inobj$Location_short)>50]))
inobj$Location_short=droplevels(inobj$Location_short)
if(is.null(GPorder)){
  GPorder=levels(inobj@meta.data[,"Location_short"])
}
groupColor=hue_pal()(length(GPorder))
names(groupColor)=GPorder

inobj=RunUMAP(inobj, dims = 1:50, reduction = "pca_harmony", reduction.name = "umap",
              n.neighbors = 100, min.dist = 0.15, spread = 1.2, verbose = FALSE)

mycolortmp=sample(sciColor_filter$hex_color, size = nCT)
show_col(mycolortmp)
names(mycolortmp)=levels(inobj@meta.data[,iCT])
reslist[[paste0(iCT,"_color")]]=mycolortmp

reslist[[paste0(iCT,"_dimplot0")]]=DimPlot_idx(inobj, group.by = iCT, cols = as.vector(mycolortmp), label.idx.size = 5, raster = F)+
  xlab("UMAP")+ylab("UMAP")+myfont

reslist[[paste0(iCT,"_roe")]]=roveRun(inobj, Group="Location_short", CT=iCT, outdir=outdir, flip = FALSE, mywidth=4.5, nCT=nCT)

reslist[[paste0(iCT,"_roeEntropy")]]=roe_entro_plot(roefile=paste0(outdir, "rove_main_age.txt"),
                   intable=table(inobj@meta.data[,iCT],inobj$Location_short), thresh=2)

if(!is.null(degDF)){
  if("cluster" %!in% colnames(degDF)){
    degDF$cluster=degDF$celltype
  }
  degDF=degDF[degDF$cluster %in% unique(inobj@meta.data[,iCT]),]
  if(!is.null(pcgeneDF)){
    degDF = degDF[degDF$gene %in% pcgeneDF$V1,]
  }
  top10 <- degDF %>% group_by(cluster) %>% top_n(6, avg_log2FC)
  subobj <- inobj
  cellID <- subsetID2(subobj@meta.data[,c(iCT,'cellid')],expected.cell = 100)
  subobj <- subset(subobj, cellid %in% cellID$Cellnames)
  subobj <- ScaleData(subobj, features = rownames(subobj))

  Idents(subobj) <- factor(Idents(subobj), levels = CTorder)
  top10$cluster <- factor(top10$cluster, levels = CTorder)
  top10 <- top10[order(top10$cluster),]
  p1 <- DoHeatmap(subobj, features = top10$gene, size = 3, group.by = iCT, group.colors = as.vector(mycolortmp),angle = 30,group.bar.height = 0.05)
  p1[[1]]$layers[[2]]$show.legend=F
  p=p1[[1]]+scale_fill_gradient2(low = "#141464", mid = "white", high = "#78050F", midpoint = 0, limits=c(-2.5,2.5))
  reslist[[paste0(iCT,"_doheatmap")]]=p
}

if(!is.null(funDF)){
  ###
  funDF=funDF[funDF$celltype %in% unique(inobj@meta.data[,iCT]),]
  df_GOs=funDF

  df_GOs_top1=df_GOs %>% group_by(celltype) %>% top_n(n = 1, wt=Factor*(-1*log10(p.adjust)))
  df_GOs_top3=df_GOs %>% group_by(celltype) %>% top_n(n = 3, wt=Factor*(-1*log10(p.adjust)))
  myterms=table(df_GOs_top3$Description)[table(df_GOs_top3$Description)>=2 & table(df_GOs_top3$Description)<=5]
  df_GOs_plotdf=df_GOs_top3[df_GOs_top3$Description %in% names(myterms),]
  df_GOs_plotdf=rbind(df_GOs_plotdf, df_GOs_top1) %>% unique()
  df_GOs_plotdf_matrix=df_GOs_plotdf %>% tidyr::pivot_wider(names_from = "Description", id_cols = "celltype", values_from = "Factor")
  rownnn=df_GOs_plotdf_matrix$celltype
  df_GOs_plotdf_matrix=df_GOs_plotdf_matrix[,2:dim(df_GOs_plotdf_matrix)[2]] %>% as.matrix()
  rownames(df_GOs_plotdf_matrix)=rownnn
  df_GOs_plotdf_matrix[is.na(df_GOs_plotdf_matrix)]=0
  dispMatrix=t(df_GOs_plotdf_matrix)
  dispMatrix=round(dispMatrix,digits = 3)
  dispMatrix[dispMatrix==0]=""
  reslist[[paste0(iCT,"_gofun_heatmap")]]=
    pheatmap(t(df_GOs_plotdf_matrix), name = "Factor", show_colnames = T, show_rownames = T,
           cluster_rows = T, cluster_cols = F,number_color ="black", display_numbers = dispMatrix, color = viridis_pal()(100), border_color = NA, number_format = "%.2f")

  ###
}

### function intergroup

}
