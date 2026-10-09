# ============================================================================
# Organ-wise macrophage and dendritic-cell subtype refinement
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Clustering and cell type annotation"
# Source     : extracted from Analysis240525.r (L4959-5457)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
OutputDir="./AnaRes20240729/"
dir.create(OutputDir, recursive = T)

### inobj (inobj=subdataList$myl_mp)
MP_resident_markers=c("FOLR2","TIMD4","LYVE1","CD74","CCR2","MRC1", "CD163", "SIGLEC1","IGF1","PDGFC","CD68","CD14")
###T cells? "CD69","RUNX3","EGR1","ITGAE","NR4A1"
inobj=fastRefine(inobj, dims = 1:16, resol = 0.3, nnei = 60, mindist =0.01, spread = 1.2)
DimPlot(inobj, group.by = "MYL_SubClass00", reduction = "umap",label = T)
DimPlot(inobj, group.by = "seurat_clusters", reduction = "umap",label = T)

inobj$MYL_SubClass01=inobj$MYL_SubClass00
inobj=modifyAno(inobj, anoCol = "MYL_SubClass01", oldID = "^Mph03_HLA$", newID = "DC01_CD1C")
inobj=modifyAno(inobj, anoCol = "MYL_SubClass01", oldID = "^Mono02_INHBA$", newID = "Mph03_INHBA")
inobj=modifyAno(inobj, anoCol = "MYL_SubClass01", oldID = "^Mono03_THBS1$", newID = "Mono02_THBS1")
MPrefine_list=AnnoRefine(inobj, refine_list=refine_list, ct="MYL_SubClass01", gp="Location_short",
                         minCellperGP=100, gp2refine=c("Skin","Lung", "Spleen", "Liver"))
### after AnnoRefine running and subdata update, MPrefine_list[[1]] was removed [SO, no MPrefine_list[[2]] exist now!]
subdata$MYL_SubClass02=subdata$MYL_SubClass01
#tissue-resident macrophages. Some of the most well-known include Langerhans cells in the skin, alveolar macrophages in the lung, Kupffer cells in the liver, microglia in the brain, and red pulp macrophages in the spleen
### Liver ###########################################################
tmp_subdata=MPrefine_list[[2]]$Liver
DotPlot(tmp_subdata, features = c("CD1C","FCN1","TIMD4", "VSIG4", "CD209", "CLEC4G","SUCNR1","MARCO","CD163","CD5L"),
        group.by = "MYL_SubClass01") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass01", reduction = "umap",label = T)
### 0,6(KC), 2(moMP)
tmp_subdata_list=list()
tmp_subdata_list[["KC01_TIMD4"]]=c(0,6)
tmp_subdata_list[["Mph02_STARD13"]]=c(2)
tmp_subdata_list[["DC01_CD1C"]]=c(7)
tmp_subdata_list[["Mono01_CD14"]]=c(1,3,4)
tmp_subdata_list[["Mono02_THBS1"]]=c(5)
###
for(x in 1:length(tmp_subdata_list)){
subdata=modifyAno(subdata, anoCol = "MYL_SubClass02",
                cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                newID = names(tmp_subdata_list[x]))
}
### Lung ###########################################################
tmp_subdata=MPrefine_list[[2]]$Lung
DotPlot(tmp_subdata, features = c("CD1C","FCN1","ITGAX", "MRC1", "SIGLEC1","MARCO","CD163","CD5L"),
        group.by = "MYL_SubClass01") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass01", reduction = "umap",label = T)
tmp_subdata_list=list()
tmp_subdata_list[["Mph01_AGRP"]]=c(0,1)
tmp_subdata_list[["AM01_ITGAX"]]=c(2)
tmp_subdata_list[["Mono01_CD14"]]=c(3)
###
for(x in 1:length(tmp_subdata_list)){
  subdata=modifyAno(subdata, anoCol = "MYL_SubClass02",
                    cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                    newID = names(tmp_subdata_list[x]))
}
### Spleen ###########################################################
tmp_subdata=MPrefine_list[[2]]$Spleen
DotPlot(tmp_subdata, features = c("CD1C","SPIC","HMOX1","FCN1","ITGAX", "MRC1", "ITGAM","MARCO","CD163","ADGRE1","CD5L"),
        group.by = "MYL_SubClass01") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass01", reduction = "umap",label = T)
### 2,4(resident RPM,red pulp MP), 0、6(moMP)
tmp_subdata_list=list()
tmp_subdata_list[["Mph04_CD5L"]]=c(0,6)
tmp_subdata_list[["RPM01_VCAM1"]]=c(2)
tmp_subdata_list[["RPM02_ITGA9"]]=c(4)
tmp_subdata_list[["Mono03_VCAN"]]=c(3,5)
tmp_subdata_list[["Mono01_CD14"]]=c(1)
tmp_subdata_list[["rmv_spleen"]]=c(7)
###
for(x in 1:length(tmp_subdata_list)){
  subdata=modifyAno(subdata, anoCol = "MYL_SubClass02",
                    cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                    newID = names(tmp_subdata_list[x]))
}

## Skin ###########################################################
tmp_subdata=MPrefine_list[[2]]$Skin
DotPlot(tmp_subdata, features =  c("CD1C","CD207","KRT5","STARD13","FCN1","ITGAX", "MRC1", "ITGAM","MARCO","CD163","CD14","CD5L"),
        group.by = "MYL_SubClass01") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass01", reduction = "umap",label = T)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mono01_CD14"]]=c(0)
tmp_subdata_list[["LC01_CD207"]]=c(1)
tmp_subdata_list[["Mph02_STARD13"]]=c(2)
tmp_subdata_list[["rmv_skin"]]=c(3)
###
for(x in 1:length(tmp_subdata_list)){
  subdata=modifyAno(subdata, anoCol = "MYL_SubClass02",
                    cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                    newID = names(tmp_subdata_list[x]))
}
## Intestine ###########################################################
tmp_subdata=MPrefine_list[[1]]$Intestine
DotPlot(tmp_subdata, features =  c("CD1C","STARD13","AGRP","FCN1","CD14","VCAN","INHBA","ITGAX", "MRC1", "ITGAM","MARCO","CD163","CD5L","TIMD4","CD4","LYVE1","FOLR2","ADAMDEC1","MEF2A"),
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
###
tmp_subdata_list=list()
tmp_subdata_list[["MphX_Intestine.C1QB"]]=c(0,1)   ###or MphX_SELENOP (M2 pheonotype related)
tmp_subdata_list[["GM01_CD4"]]=c(3,5,8)
tmp_subdata_list[["Mono01_CD14"]]=c(2,6)
tmp_subdata_list[["DC01_CD1C"]]=c(9)
tmp_subdata_list[["MphX_CCL3"]]=c(4)  ### monocyte derive MP
tmp_subdata_list[["rmv_Intestine"]]=c(7,10,11)
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                    cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                    newID = names(tmp_subdata_list[x]))
}

### Next organ?
## Airway? UDT? BS? Artery? ... ###########################################################
## UDT/Airway/Artery/BS/UB/Kidney/Ovary.Testis/VD.FU
inobj=subdataList$myl_mp
MPrefine_list_add=list()
gp2refine=c("UDT","Airway","Artery","BS","UB","Kidney","Ovary.Testis","VD.FU")
for(x in gp2refine){
  print(paste0(x,"..."))
  xv=strsplit(x, split = "\\.") %>% unlist()
  MPrefine_list_add[[x]]=fastRefine(inobj[,inobj@meta.data[,"Organ"] %in% xv],
                                    dims = 1:12, resol = 0.2, nnei = 50, mindist =0.01, spread = 0.5)
}
MP_selMarkers=c("CD1C","STARD13","AGRP","INHBA","SIGLEC8","SELENOP",
                "FCN1","CD14","VCAN","THBS1","LYPD2",
                "ITGAX", "MRC1", "ITGAM","MARCO","C1QA","CD163","ADGRE1","SIGLEC1","SPIC","CD5L","TIMD4","CD4","LYVE1","FOLR2","ADAMDEC1","MEF2A","ITGA9","VCAM1","CD207")
## UDT ###########################################################
tmp_subdata=MPrefine_list_add$UDT
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["MphX_UDT.C1QB"]]=c(0)   ###
tmp_subdata_list[["Mph02_STARD13"]]=c(1)
tmp_subdata_list[["Mono01_CD14"]]=c(2,3)
tmp_subdata_list[["DC01_CD1C"]]=c(4)
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                    cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                    newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_UDT.C1QB"), ident.2 = c("Mph01_AGRP"), group.by = "MYL_SubClass04",
                test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## Airway ###########################################################
tmp_subdata=MPrefine_list_add$Airway
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(0,6)   ###
tmp_subdata_list[["MphX_Airway.C1QB"]]=c(1,3)   ### almost No DEGs...
tmp_subdata_list[["Mono01_CD14"]]=c(2,5)
tmp_subdata_list[["MphX_Airway.LPL"]]=c(4)  ### rename as LAM (lipid related)
tmp_subdata_list[["DC01_CD1C"]]=c(7)  ### high in CD207 ?
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_Airway.C1QB"), ident.2 = c("Mph01_AGRP"),
                group.by = "MYL_SubClass04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## Artery ###########################################################
tmp_subdata=MPrefine_list_add$Artery
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(1)   ###
tmp_subdata_list[["MphX_Artery.C1QB"]]=c(2)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(0,3)
tmp_subdata_list[["DC01_CD1C"]]=c(4)  ###
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_Artery.C1QB"), #ident.2 = c("Mph01_AGRP"),
                group.by = "MYL_SubClass04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## BS ###########################################################
tmp_subdata=MPrefine_list_add$BS
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(0,5,6)   ###
tmp_subdata_list[["MphX_BS.C1QB"]]=c(1)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(2,3)
tmp_subdata_list[["DC01_CD1C"]]=c(4)  ###
tmp_subdata_list[["rmv_BS"]]=c(7)
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_BS.C1QB"), #ident.2 = c("Mph01_AGRP"),
                group.by = "MYL_SubClass04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## UB ###########################################################
tmp_subdata=MPrefine_list_add$UB
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(0)   ###
tmp_subdata_list[["MphX_UB.C1QB"]]=c(1)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(2,3)
tmp_subdata_list[["DC01_CD1C"]]=c(4)  ###
tmp_subdata_list[["Mph05_SIGLEC8"]]=c(5)  ###! Overalp good; but SIGLEC8 exp is LOW!
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_UB.C1QB"), #ident.2 = c("Mph01_AGRP"),
                group.by = "MYL_SubClass04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## Kidney ###########################################################
tmp_subdata=MPrefine_list_add$Kidney
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(0,2)   ###
tmp_subdata_list[["DC02_CEACAM4"]]=c(3)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(1)
tmp_subdata_list[["DC01_CD1C"]]=c(4)  ###
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
xxx=FindMarkers(inobj,ident.1 = c("MphX_Kidney.C1QB"), #ident.2 = c("Mph01_AGRP"),
                group.by = "MYL_SubClass04", test.use = "t", max.cells.per.ident = 2000, only.pos = TRUE, min.pct = 0.25)
## Ovary.Testis ###########################################################
tmp_subdata=MPrefine_list_add$Ovary.Testis
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(0)   ###
tmp_subdata_list[["MphX_OvaryTestis.C1QB"]]=c(1)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(2)
tmp_subdata_list[["MonoX_OvaryTestis"]]=c(3)  ###
###
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid,
                  newID = names(tmp_subdata_list[x]))
}
#xxx=FindMarkers(inobj,ident.1 = c("MphX_Kidney.C1QB"), #ident.2 = c("Mph01_AGRP"),
## VD.FU ###########################################################
tmp_subdata=MPrefine_list_add$VD.FU
tmp_subdata$MYL_SubClass03=droplevels(tmp_subdata$MYL_SubClass03)
tmp_subdata$Location_short=droplevels(tmp_subdata$Location_short)
DotPlot(tmp_subdata, features =  MP_selMarkers,
        group.by = "MYL_SubClass03") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+xlab("")+ylab("")
DimPlot(tmp_subdata, group.by = "MYL_SubClass03", reduction = "umap",label = T)
table(tmp_subdata$MYL_SubClass03,tmp_subdata$Location_short)
###
tmp_subdata_list=list()
tmp_subdata_list[["Mph02_STARD13"]]=c(2)   ###
tmp_subdata_list[["MphX_VDFU.C1QB"]]=c(0)   ###
tmp_subdata_list[["Mono01_CD14"]]=c(1)
###
negVec=subset(tmp_subdata, MYL_SubClass03 %in% c("DC01_CD1C"))$cellid %>% as.vector()
for(x in 1:length(tmp_subdata_list)){
  inobj=modifyAno(inobj, anoCol = "MYL_SubClass04",
                  cellid = setdiff(subset(tmp_subdata, seurat_clusters %in% tmp_subdata_list[[x]])$cellid, negVec),
                  newID = names(tmp_subdata_list[x]))
}
#xxx=FindMarkers(inobj,ident.1 = c("MphX_Kidney.C1QB"), #ident.2 = c("Mph01_AGRP"),

inobj_sub=subset(inobj, MYL_SubClass04_bak %!in% c("rmv_BS","rmv_Intestine","rmv_skin","rmv_spleen","rmv_UDT"))
inobj_sub$MYL_SubClass04_bak=droplevels(inobj_sub$MYL_SubClass04_bak)
xxx=colPercents(table(inobj_sub$MYL_SubClass04_bak, inobj_sub$Location_short, inobj_sub$Donor)) %>% as.data.frame() %>% head(25)
write.xlsx(xxx %>% as.data.frame(), file = paste0("./MPsub_byDonorLocation_allcombined.xlsx"),colNames = TRUE, rowNames = TRUE)

##################################### NOT done yet in 0825: subdataList$myl_mp=inobj
###!!! added DONE in 20241011 subdataList$myl_mp=inobj (with "MYL_SubClass04", MYL_SubClass04_bak,ClassActMP.01, AlterActMP.01 ) have DONE!!!!
