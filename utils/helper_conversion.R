# ============================================================================
# Data conversion and small text helpers (h5ad2seu, CorTable, toStr, text2vector)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Source     : extracted from Analysis240525.r (L21-60, L1593-1637, L6797-6803)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
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

CorTable=function(dfin, grp="group", cluster="cluster", gene="gene", JC=TRUE, ...){

grp_item=unique(dfin[,grp])
df1=dfin[dfin[,grp]==grp_item[1],]
df2=dfin[dfin[,grp]==grp_item[2],]

clu1=as.vector(unique(df1[,cluster]))
clu2=as.vector(unique(df2[,cluster]))

out.matrix <- matrix(0, nrow = length(clu1), ncol = length(clu2))

for(i in 1:length(clu1)){
  for(j in 1:length(clu2)){
    #i=1
    #j=1
    geneS1=as.vector(unique(df1[df1[,cluster]==clu1[i],gene]))
    geneS2=as.vector(unique(df2[df2[,cluster]==clu2[j],gene]))
    if(JC){
      jacDis=round(length(intersect(geneS1,geneS2))/length(union(geneS1,geneS2)),digits = 3)
    }else{
      jacDis=length(intersect(geneS1,geneS2))
    }
    out.matrix[i,j]=jacDis
  }
}

rownames(out.matrix)=clu1
colnames(out.matrix)=clu2

return(out.matrix)

}

toStr=function(df,...){
tc <- textConnection("str", "w")
sink(tc)
print(df)
sink()
close(tc)

str=gsub(pattern = "\\s+.+", replacement = "", str, perl = T)
str <- paste0(" ",str," ",collapse="\n")
cat(str)
}

text2vector=function(intext, sep="\n"){
  xxx=unlist(strsplit(intext, split=sep, perl = TRUE))
  xxx=gsub("^\\s+","",xxx)
  xxx=gsub("\\s+$","",xxx)
  xxx=xxx[nchar(xxx)!=0]
  return(xxx)
}
