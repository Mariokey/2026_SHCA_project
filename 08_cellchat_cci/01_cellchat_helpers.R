# ============================================================================
# CellChat wrappers (runCCI_perData, matrix_expand, pickDotCCI)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Cell-cell interactions"
# Source     : extracted from Analysis240525.r (L7417-7568, L7569-7719)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
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
runCCI_perData=function(obj, spatialRun=FALSE, inlist=NULL,
                      group="NC_LV", groupby=NULL, indir="DATA_DIR", drawOnly=FALSE){

###group: levels in groupby or prefix mark
runCCI=function(cellchat_seurat, spatial=FALSE, inST=NULL){

  inmatrix=cellchat_seurat@assays$RNA$data

  if(!spatial){
    cellchat <- createCellChat(object = inmatrix,
                               meta = cellchat_seurat@meta.data,
                               group.by = "labels")   ###cellchat_seurat@assays$RNA@data
  }else{
    spatial.locs=inST[["loc"]]
    spatial.factors=inST[["factor"]]
    cellchat <- createCellChat(object = inmatrix, meta = cellchat_seurat@meta.data, group.by = "labels",
                               datatype = "spatial", coordinates = spatial.locs,
                               spatial.factors = spatial.factors)
  }
  cellchat@DB <- CellChatDB   ### CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if(!spatial){
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)    ###type = "truncatedMean" and trim = 0.1
  }else{
    d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
    min_d=min(d.spatial[d.spatial!=0])
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean",
                                  trim = 0, k.min = 5, nboot = 20,
                                  distance.use = TRUE, interaction.range = 200, scale.distance = 1.5/min_d,
                                  contact.dependent = FALSE, contact.range = NULL)
  }
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

dir.create(indir, recursive = T)
rds_file=paste0(indir,"/",paste0(group,collapse = ""),".CCI.rds")

if(!spatialRun){inlist=NULL}

if(!drawOnly){
  if(!is.null(groupby)){
    tVobj=obj[, obj@meta.data[,groupby] %in% c(group)] #subset(obj, Sample %in% c(group))  #### Sample: "RO"  "R3d" "SO"  "S3d"
    tVobj@meta.data[,groupby]=droplevels(tVobj@meta.data[,groupby])
  }else{
    tVobj=obj
  }
  tVobj$labels=tVobj$CCItype
  tVobj$labels = droplevels(tVobj$labels, exclude = setdiff(levels(tVobj$labels),unique(tVobj$labels)))

  if(!is.null(inlist)){
    if(inlist[["mode"]]=="stereo"){
      spatial.locs=tVobj@reductions$spatial@cell.embeddings %>% as.data.frame()
      colnames(spatial.locs)=c("imagecol","imagerow")
    }else{
      sampleDir_base=names(tVobj@images)
      for(x in sampleDir_base[sampleDir_base %!in% group]){ tVobj@images[[x]]=NULL }
      spatial.locs = Seurat::GetTissueCoordinates(tVobj, scale = NULL, cols = c("imagerow", "imagecol"))
    }
    inlist[['loc']]=spatial.locs
  }

  tVobj=runCCI(tVobj, spatial=spatialRun, inST=inlist)
  saveRDS(tVobj, file = rds_file)
}else{
  tVobj=readRDS(rds_file)
}

if(max(tVobj@net$count)>0){

  group=paste0(group,collapse = "")

  groupSize <- as.numeric(table(tVobj@idents))
  tdf=as.data.frame(tVobj@net$count)
  A=pheatmap(as.matrix(tdf),cluster_rows = F,cluster_cols = F, display_numbers = as.matrix(tdf))
  A=Figgo(text = paste0(group,".","Heatmap_Ninter",".CCI"), dir_img = indir, obj = A, w = 7, h = 5)

  CairoPDF(paste0(indir,"/",group,".","Circle_Ninter",".CCI",'.pdf'), 6, 6.5)
  netVisual_circle(tVobj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  CairoPDF(paste0(indir,"/",group,".","Circle_Nweight",".CCI",'.pdf'), 6, 6.5)
  netVisual_circle(tVobj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  if(spatialRun){
    CairoPDF(paste0(indir,"/",group,".","Spatial_Nweight",".CCI",'.pdf'), 7, 5.5)
    xxx=netVisual_spatial(tVobj@net$weight, tVobj@images$coordinates, tVobj@meta,
                          vertex.weight = groupSize, weight.scale = T, sample.use = group,
                          point.size = 0.3,edge.width.max = 2, vertex.weight.max = median(groupSize),
                          vertex.size.max = 8, alpha.image = 0.2, vertex.label.cex = 0)
    plot(xxx)
    dev.off()
  }
}
}

matrix_expand=function(xxxh, identV, rename_list=NULL){

xxxm=matrix(0.0000001,nrow = length(identV), ncol = length(identV),dimnames = list(identV,identV))
if(!is.null(rename_list)){
  for(i in 1:length(rename_list)){
    rownames(xxxh)=gsub(paste0('^',rename_list[[i]][1],'$'),rename_list[[i]][2],rownames(xxxh))
    colnames(xxxh)=gsub(paste0('^',rename_list[[i]][1],'$'),rename_list[[i]][2],colnames(xxxh))
  }
}
for(x in identV){
  for(y in identV){
    xxxm[x,y]=ifelse(x %in% rownames(xxxh) & y %in% rownames(xxxh),xxxh[x,y],0.0000001)
  }
}

return(xxxm)
}
pickDotCCI=function(t12merge, fCT, LRs=NULL, Source=TRUE, Target=TRUE, selfrmv=TRUE, flip=TRUE, ...){

idx=levels(t12merge@idents$joint)
bulb=netVisual_bubble(t12merge, sources.use = c(1:length(idx)), targets.use = c(1:length(idx)),
                      comparison = c(1,2), angle.x = 45, return.data=T)
bulbtmp=bulb$gg.obj
intergroups=unique(bulbtmp$data$dataset) %>% as.vector()

myLRs=levels(bulbtmp$data$interaction_name)
myLRs=myLRs[!is.na(myLRs)]
mylevel=apply(expand.grid(unique(bulbtmp$data$group.names),
                          paste(" (",unique(bulbtmp$data$dataset),")", sep = "")), 1, paste, collapse="")
zzz=mylevel
mycolor=ifelse(grepl("3D",zzz), "#F8766D", ifelse(grepl("4W",zzz), "Orange", ifelse(grepl("B",zzz), "#619CFF", "#00BFC4")))
names(mycolor)=zzz

bulbtmp$data$interaction_name=as.character(bulbtmp$data$interaction_name)
bulbtmp$data=bulbtmp$data[!is.na(bulbtmp$data$interaction_name),]
bulbtmp$data$source.target=as.character(bulbtmp$data$source.target)
bulbtmp$data=bulbtmp$data[!is.na(bulbtmp$data$source.target),]

if(!is.null(LRs)){
  LRpool=c()
  for(v in LRs){
    LRpool=c(LRpool, unique(bulbtmp$data$interaction_name[
      grepl(v, bulbtmp$data$interaction_name)])) ### %>% as.vector()))
  }
  bulbtmp$data=bulbtmp$data[bulbtmp$data$interaction_name %in% LRpool,]
}

#                           (bulbtmp$data$source %!in% c(AllNeu) & bulbtmp$data$target %in% c(fCT)),]
bulbtmp$data=bulbtmp$data[bulbtmp$data$source %in% c(fCT) | bulbtmp$data$target %in% c(fCT), ]
if(!Source){
  if(selfrmv){
    bulbtmp$data=bulbtmp$data[bulbtmp$data$source %!in% c(fCT),]
  }else{
    bulbtmp$data=bulbtmp$data[(bulbtmp$data$target %in% c(fCT)),]
  }
}
if(!Target){
  if(selfrmv){
    bulbtmp$data=bulbtmp$data[bulbtmp$data$target %!in% c(fCT),]
  }else{
    bulbtmp$data=bulbtmp$data[(bulbtmp$data$source %in% c(fCT)),]
  }
}

yyy_ori=table(bulbtmp$data$group.name, bulbtmp$data$dataset) #%>% as.data.frame()
yyy=matrix(0,nrow = nrow(yyy_ori), ncol = length(intergroups),dimnames = list(rownames(yyy_ori),intergroups))
for(x in rownames(yyy_ori)){
  for(y in intergroups){
    yyy[x,y]=ifelse(x %in% rownames(yyy_ori) & y %in% colnames(yyy_ori),yyy_ori[x,y],0)
  }
}
yyy=yyy[MatrixGenerics::rowSums(yyy)>0,]
yyy_inter_level=apply(expand.grid(rownames(yyy),
                                  paste(" (",colnames(yyy),")", sep = "")), 1, paste, collapse="")
yyy_inter_level=sort(yyy_inter_level)
keepLab=rownames(yyy)
bulbtmp$data=bulbtmp$data[bulbtmp$data$group.name %in% keepLab,]
add_level=setdiff(yyy_inter_level, unique(bulbtmp$data$source.target))
if(FALSE & length(add_level)>0){
  add_level=add_level[1]
  add_level_v= add_level %>% strsplit(split = " \\(|\\)", perl = TRUE) %>% unlist()
  bulbtmp$data[nrow(bulbtmp$data)+1,]=NA #bulbtmp$data[nrow(bulbtmp$data),]
  bulbtmp$data[nrow(bulbtmp$data),"source.target"]=add_level
  bulbtmp$data[nrow(bulbtmp$data),"interaction_name_2"]=as.character(bulbtmp$data$interaction_name_2[1])
  bulbtmp$data[nrow(bulbtmp$data),"group.names"]=add_level_v[1]
  bulbtmp$data[nrow(bulbtmp$data),"dataset"]=add_level_v[2]
  bulbtmp$data[nrow(bulbtmp$data),"prob"]=1
  bulbtmp$data[nrow(bulbtmp$data),"pval"]=0
}

#
#                  keepLab]
mylevel=yyy_inter_level
bulbtmp$data$source.target=factor(bulbtmp$data$source.target, levels=mylevel)
myLRs=myLRs[myLRs %in% unique(bulbtmp$data$interaction_name)]
bulbtmp$data$interaction_name=factor(bulbtmp$data$interaction_name, levels = myLRs)

if(length(unique(bulbtmp$data$pval))==1){
  #bulbtmp$data[1,]$pval=if_else(xxxx<=2, xxxx+1, xxxx-1)
  xxxx=bulbtmp$data[bulbtmp$data$prob.original==min(bulbtmp$data$prob.original),]$pval
  bulbtmp$data[bulbtmp$data$prob.original==min(bulbtmp$data$prob.original),]$pval=xxxx-1
}

bulbtmp$theme$axis.text.y$colour=NULL
bulbtmp$theme$axis.text.x$colour=NULL
bulbtmp$layers[[3]]=NULL
bulbtmp$layers[[2]]=NULL

if(flip){
  bulbtmp=bulbtmp+coord_flip()
  bulbtmp$theme$axis.text.y$colour=mycolor[levels(bulbtmp$data$source.target)]
  bulbtmp=bulbtmp+scale_x_discrete(drop = FALSE)
}else{
  bulbtmp$theme$axis.text.x$colour=mycolor[levels(bulbtmp$data$source.target)]
  bulbtmp=bulbtmp+scale_x_discrete(drop = FALSE)
}

print(head(bulbtmp$data))
return(bulbtmp)
}
