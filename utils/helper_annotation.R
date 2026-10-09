# ============================================================================
# Annotation helpers (annotate_cell_types, fastRefine, AnnoRefine, addNewMetaColumn)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Source     : extracted from Analysis240525.r (L1294-1309, L1639-1646, L6857-6889, L7417-7451)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
annotate_cell_types <- function(seurat_obj, cell_type_vector, newColID="NKSubtypes01") {

# Check if all active.ident in the vector are present in the Seurat object
#}

# Add the cell type annotations to the Seurat object's metadata
seurat_obj@meta.data[,newColID] <- cell_type_vector[seurat_obj@active.ident]
seurat_obj@meta.data[,newColID]=factor(seurat_obj@meta.data[,newColID], levels = sort(unique(as.vector(cell_type_vector))))

return(seurat_obj)
}

fastRefine=function(inobj, reduct = "pca_harmony", dims = 1:15, resol = 1.2, alg = 2, nnei = 120, mindist = 0.15, spread = 1.5){
### inobj already with "pca_harmony"
inobj <- FindNeighbors(inobj, reduction = reduct, dims = dims, k.param = nnei, verbose = FALSE)
inobj <- RunUMAP(inobj, dims = dims, reduction=reduct, reduction.name = "umap",
                 n.neighbors = nnei, min.dist = mindist, spread = spread, verbose = FALSE)
inobj <- FindClusters(inobj, resolution = resol, algorithm = alg, verbose = FALSE)
return(inobj)
}

AnnoRefine=function(obj, refine_list=NULL, ct="MYL_SubClass00", gp="Location_short", minCellperGP=100, gp2refine=NULL){

  tmpdatalist=list()

  obj$Organ=obj@meta.data[,gp]
  for(i in 1:length(refine_list)){
    for(j in 1:length(refine_list[[i]])){
    obj=modifyAno(obj, anoCol = "Organ", oldID = paste0("^",refine_list[[i]][j],"$"), newID = names(refine_list[i]))
    }
  }
  obj$Organ=factor(obj$Organ, levels = rev(names(refine_list)))

  organ_table=table(obj$Organ)
  organ_table=organ_table[organ_table>=minCellperGP]
  print(organ_table)
  obj=obj[,obj@meta.data[,"Organ"] %in% names(organ_table)]

  if(is.null(gp2refine)){
    gp2refine=names(organ_table)
  }else{
    gp2refine=intersect(names(organ_table),gp2refine)
  }
  for(x in gp2refine){
    print(paste0(x,"..."))
  tmpdatalist[[x]]=fastRefine(obj[,obj@meta.data[,"Organ"] %in% x],
                              dims = 1:12, resol = 0.1, nnei = 50, mindist =0.01, spread = 1)
  }

  return(list(obj, tmpdatalist))
}

addNewMetaColumn=function(Mobj, tmplist, nameV=NULL, forceAdd=FALSE, forceAddname=NULL, ...){

if(forceAdd){
  Mobj$OKtmpColumnXXX=Mobj@meta.data[,forceAddname]
  newName=forceAddname
}else{
  Mobj$OKtmpColumnXXX=NA
  newName="OKtmpColumnXXX"
  if(!is.null(nameV)){
    newName=nameV
  }
}

for(i in 1:length(tmplist)){
  tmp_obj=names(tmplist[i])
  tmp_colnames=tmplist[i][[1]]
  subaddobj=eval(as.name(tmp_obj))
  for(j in 1:length(tmp_colnames)){
    tmpAnno2=Mobj$OKtmpColumnXXX
    tmpAnno3=as.character(subaddobj@meta.data[,tmp_colnames[j]])
    names(tmpAnno3)=colnames(subaddobj)

    tmpAnno2[colnames(subaddobj)]=tmpAnno3
    Mobj=AddMetaData(Mobj, metadata = tmpAnno2, col.name = newName)
    Mobj$OKtmpColumnXXX=Mobj@meta.data[,newName]
  }
}

if(newName!="OKtmpColumnXXX"){
  Mobj$OKtmpColumnXXX=NULL
}
return(Mobj)

}
