# ============================================================================
# tradeSeq gene expression patterns along the slingshot pseudotime (tradeS)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Trajectory analysis"
# Source     : extracted from Analysis240525.r (L4467-4639)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
tradeS=function(sce, samplesize=0.05, SvE=FALSE,  Asso=TRUE, outdirprefix="DATA_DIR", ...){

  ### Tradeseq analysis : within 30min < 4000 cells
  library(tradeSeq)
  library(RColorBrewer)
  library(SingleCellExperiment)
  library(slingshot)
  library(scran)
  library(org.Hs.eg.db)

  tradeResList=list()

BPPARAM = BiocParallel::bpparam()
BPPARAM$workers=6

sce1k=sce[, sample(colnames(sce), size=round(samplesize*length(colnames(sce))))]

tseq_counts=sce1k@assays@data@listData$counts
tseq_pseudotime <- slingPseudotime(sce, na = FALSE)  ####!!! na=FALSE: arclength along each curve will be returned for NA cells.
if(dim(tseq_pseudotime)[2]==1){
  tseq_pseudotime=data.frame(Lineage1=as.vector(tseq_pseudotime[colnames(sce1k),]),
                             row.names = names(tseq_pseudotime[colnames(sce1k),]))
}else{
  tseq_pseudotime=tseq_pseudotime[colnames(sce1k),]
}
tseq_cellWeights <- slingCurveWeights(sce1k)

tradeResList[["tradeS_sce1k"]]=sce1k
tradeResList[["tradeS_sce1k_time"]]=tseq_pseudotime

set.seed(7)
x <- org.Hs.egGENETYPE
mapped_type <- mappedkeys(x)
xx <- as.list(x[mapped_type])
pcgene_id=names(xx[xx=="protein-coding"])
y <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(y)
yy <- as.list(y[mapped_genes])
pcgene_symbol=as.vector(unlist(yy[pcgene_id]))

variable_genes <-
  sce1k %>%
  scran::modelGeneVar() %>%
  scran::getTopHVGs(prop=0.1)

tseq_selectGenes=variable_genes  ### or from obj seurat data
tseq_selectGenes=tseq_selectGenes[tseq_selectGenes %in% pcgene_symbol]
tseq_sce <- fitGAM(counts = tseq_counts, pseudotime = tseq_pseudotime, cellWeights = tseq_cellWeights, #conditions = factor(obj$Group),
                   parallel = TRUE, BPPARAM = BPPARAM, control = mgcv::gam.control(),
                   genes=tseq_selectGenes, nknots = 7, verbose = T)
tradeResList[["tradeS_sce1k_tseq_sce"]]=tseq_sce

num_Vgene_used=length(tseq_selectGenes)
print(paste0("Number of input variable genes for tradeS: ", num_Vgene_used))
table(rowData(tseq_sce)$tradeSeq$converged)

if(Asso){
  assoRes <- associationTest(tseq_sce, lineages=TRUE)
  tradeResList[["tradeS_assoRes"]]=assoRes
}

if(SvE){
  startRes <- startVsEndTest(tseq_sce, lineages=TRUE)  ##pseudotimeValues = c(0, 1)  ##for large data l2fc = log2(2) to reduce FalsePositiveRate
  tradeResList[["tradeS_startRes"]]=startRes
}

  if(Asso){
    tempAll=as.data.frame(t(tseq_counts[rownames(assoRes),rownames(tseq_pseudotime)]))
    tempAll$cellid=rownames(tempAll)

    for(mark in 1:dim(tseq_pseudotime)[2]){

      targetL=assoRes[assoRes[,paste0("pvalue_",mark)]<0.05,]
      temp=data.frame(cellid=rownames(tseq_pseudotime), #names(selectTime),
                      Lineage=as.vector(tseq_pseudotime[,paste0("Lineage",mark)]))  #as.vector(selectTime))
      intersctGenes=intersect(colnames(tempAll),rownames(targetL))
      temp=temp[order(temp[,"Lineage"], decreasing = F),]
      temp=left_join(temp, tempAll[,c("cellid", intersctGenes)], by="cellid")

      tradeResList[[paste0("L",mark,"_rawExp")]]=temp

      nana=data.frame(Lineage=temp$Lineage)
      for(i in 3:dim(temp)[2]){
        temp$gene=temp[,i]
        loessmod <- loess(gene ~ Lineage, data = temp, span = 1)
        smoothed_data <- predict(loessmod)
        smooth_scale <- rescale(smoothed_data, to=c(0,100))
        nana <- cbind(nana, smooth_scale)
      }
      temp$gene=NULL
      temp$cellid=NULL
      colnames(nana)=colnames(temp)
      #tradeResList[[paste0("L",mark,"_smoothExp")]]=nana
tmpheat=t(as.matrix(nana[,-1]))
colnames(tmpheat)=nana$Lineage
gene_weights=c()
pseudotime=as.numeric(colnames(tmpheat))
for(i in 1:dim(tmpheat)[1]){
  gene_val=(tmpheat[i,]+1)/sum(tmpheat[i,])
  gene_weights=c(gene_weights, sum(gene_val*pseudotime))
}
names(gene_weights)=rownames(tmpheat)
gene_weights=sort(gene_weights,decreasing = T)
AAA=pheatmap(tmpheat[names(gene_weights),], color = viridisLite::viridis(256, option = "B"), name = "Relative\nexpression",
             heatmap_legend_param = list(direction = "horizontal", title_position = "leftcenter"),
             cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
AAA=draw(AAA, heatmap_legend_side = "bottom")

tradeResList[[paste0("L",mark,"_smoothExp")]]=tmpheat[names(gene_weights),]
tradeResList[[paste0("L",mark,"_smoothExp_heatplot")]]=AAA

      original_gene_list=gene_weights   ###avg_log2FC MUST decreasing=TRUE
      original_gene_list=rescale(original_gene_list, to=c(-1,1))
      res.GSE=gseFun(original_gene_list)
      res.GSE$gse_fdf=res.GSE$gse_fdf[order(res.GSE$gse_fdf$NES, decreasing = T),]  ### !!!! IMPORTANT, NES>0 are enrichment of genes at late pseudotime
      head(res.GSE$gse_fdf)
      tradeResList[[paste0("L",mark,"_GSEres")]]=res.GSE
    }
  }

  saveRDS(tradeResList, file = paste0(outdirprefix,".tradeseqResult.rds"))
  return(tradeResList)
}
