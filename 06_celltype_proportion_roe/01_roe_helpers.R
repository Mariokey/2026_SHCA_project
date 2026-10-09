# ============================================================================
# Ro/e computation helpers (roveRun, entropy_cal, roe_entro_plot, ROIE, divMatrix)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Enrichment analysis of cell type proportion"
# Source     : extracted from Analysis240525.r (L1236-1263, L6805-6855, L14749-14794)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
roveRun=function(obj, Group="Age", CT="CT89", gorder=NULL, ctorder=NULL, flip=TRUE, outdir="DATA_DIR", mywidth=3.5, nCT=23, ...){

# Ro/e
meta <- obj@meta.data[,c(Group,CT)]
# Age
rovemat <- scPioneer::rove(x=meta[,CT],
                y=meta[,Group], display_numbers = T,
                cluster_rows = F, cluster_cols = F, plot=F)
write.table(rovemat, sep = '\t', quote = F, file = paste0(outdir,'rove_main_age.txt'))
color <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[c(1,3,5,7,9)]
lebreaks <- c(0,1,1.5,2,2.5)
graphics.off()
pdf(paste0(outdir, 'Roe_',Group,'.pdf'), mywidth, 7*(nCT/23))
if(!is.null(ctorder)){
  rovemat=rovemat[ctorder,]
}
if(!is.null(gorder)){
  rovemat=rovemat[,gorder]
}
xxx=scPioneer::plotRO(rovemat, rowname.order =ctorder, colname.order = gorder,
           cluster_cols = FALSE, cluster_rows = FALSE,
           breaks.hp = lebreaks, legend.breaks = lebreaks,
           paletteLength = 5,cellheight = 18,
           color.ramp = color, coord.flip = flip)
plot_grid(xxx$gtable, scale = 0.95)
dev.off()
return(plot_grid(xxx$gtable, scale = 0.95))
}

entropy_cal=function(zzz, rowcol="c"){

  entropy_res=c()
  if(rowcol=="r"){
    zzzr=rowPercents(zzz)[,1:dim(zzz)[2]]
    myzzz=t(zzzr)
  }
  if(rowcol=="c"){
    zzzc=colPercents(zzz)[1:dim(zzz)[1],]
    myzzz=zzzc
  }
  for(x in 1:dim(myzzz)[2]){
    entropy_res=c(entropy_res,
                  entropy.empirical(rep(100/dim(myzzz)[1],dim(myzzz)[1]), unit="log2")-
                    entropy.empirical(myzzz[,x], unit="log2"))
  }
  names(entropy_res)=colnames(myzzz)

  return(entropy_res)
}

roe_entro_plot=function(roefile="./AnaRes20240718/rove_main_age.txt", intable, thresh=2){
  tmp_entropy=entropy_cal(intable, rowcol = "c")
  yyy2=data.frame(group=names(tmp_entropy), entropy=as.vector(tmp_entropy))
  yyy2$group=factor(yyy2$group, levels = names(tmp_entropy))

  yyy1=read.table(roefile)
  yyy1$celltype=rownames(yyy1)
  yyy1=yyy1 %>% tidyr::pivot_longer(cols = 1:(dim(yyy1)[2]-1))
  colnames(yyy1)=c("celltype","group","roe")
  yyy1$group=factor(yyy1$group, levels = names(tmp_entropy))
  yyy1$text=yyy1$celltype
  yyy1[log2(yyy1$roe+1)<thresh,]$text=""

  yplot1=ggplot()+geom_bar(data=yyy2, aes(group,entropy), stat = "identity", fill="darkgrey", color="white")+
    theme_minimal_hgrid()+ylab("Entropy")+theme(axis.text.x = element_blank(),axis.title.x = element_blank())

  yplot2=ggplot()+
    geom_point(data = yyy1[yyy1$text=="",], aes(group, log2(roe+1)), color="lightgrey",
               stat = "identity", position = position_jitter(width = 0.1))+
    geom_point(data = yyy1[yyy1$text!="",], aes(group, log2(roe+1)), color="red", size=2,
               stat = "identity", position = position_jitter(width = 0.1))+
    ggrepel::geom_text_repel(data=yyy1[yyy1$text!="",], mapping=aes(group, log2(roe+1),
                                                                    label=text))+
    theme_cowplot()+ylab('Ro/e (log scale)')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),axis.title.x = element_blank())
  rplot=plot_grid(yplot1, yplot2, ncol = 1, rel_heights = c(1,3), align = "v", axis = "lrtb")

  return(rplot)
}

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1.
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}
