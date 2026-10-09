# ============================================================================
# Cell-cell communication at subtype resolution for the macrophage / CNTN4+ axis
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Cell-cell interactions"
# Source     : extracted from Analysis240525.r (L13984-14115)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
### CCI 0805

CCIdir_img_root="./CCI_sub0805/"
dir.create(CCIdir_img_root)

for(cmt in unique(LOC_weight$CMT)){

  CCIdir_img=paste0("./CCI_sub0805/",cmt,"/")
  dir.create(CCIdir_img)

  sample_oi=LOC_weight[LOC_weight$CMT %in% cmt,]$sample
  tmp_score=CT_weight[,cmt]
  tmp_score=tmp_score[tmp_score>=0.015]
  ct_oi=names(tmp_score)

rmv_ct="Prolif."
Mobj_tmp=Mobj_keep #subset(Mobj_keep, Organ %!in% c("Pancreas","Testis","VD","Ovary","FU"))
Mobj_tmp=subset(Mobj_tmp, subtypes015 %!in% c(rmv_ct))
Mobj_tmp=subset(Mobj_tmp, subtypes015 %in% c(ct_oi))
Mobj_tmp=subset(Mobj_tmp, sample %in% c(sample_oi))

Mobj_tmp$Organ=droplevels(Mobj_tmp$Organ)
Mobj_tmp$subtypes015=droplevels(Mobj_tmp$subtypes015)

Mobj_tmp$CCItype=Mobj_tmp$subtypes015
runCCI_perData(Mobj_tmp,
               spatialRun = FALSE, #inlist = st_preCfg,
               group = levels(Mobj_tmp$Organ), groupby = "Organ",
               indir = CCIdir_img, drawOnly=FALSE)

}

CCIrds=list.files(CCIdir_img)
CCIrds=CCIrds[grepl(".rds",CCIrds)]
xi=CCIrds[1]

allCCIlr_DF=data.frame()
for(xi in CCIrds){
  xi_data=readRDS(paste0(CCIdir_img,"/",xi))
  xi=gsub(".CCI.rds","",xi)

  idx=levels(xi_data@idents)
  bulb=netVisual_bubble(xi_data, sources.use = c(1:length(idx)), targets.use = c(1:length(idx)),
                        angle.x = 45, return.data=T)
  tmpDF=bulb$communication
  tmpDF$tissue=xi

  allCCIlr_DF=rbind(allCCIlr_DF, tmpDF)
}
write.xlsx(allCCIlr_DF, file = paste0(CCIdir_img,"/","All_subOI_CCIlr_DF.across_tissue.xlsx"))

tPR03=xi_data
tPR03=netAnalysis_computeCentrality(tPR03)
CellChat::netAnalysis_signalingRole_scatter(tPR03)
ggsave(paste0(CCIdir_img, "/", "netAnalysis_signalingRole_scatter.pdf"), width = 8, height = 8)

CToi="Mph03_FOLR2"
rankNet(tPR03,mode = "single",sources.use = CToi, stacked = T, do.stat = TRUE)+
  rankNet(tPR03,mode = "single",targets.use = CToi, stacked = T, do.stat = TRUE)
ggsave(paste0(CCIdir_img, "/", CToi, "_rankNet.pdf"), width = 12, height = 6)

netVisual_bubble(tPR03, sources.use = CToi)+coord_flip()
ggsave(paste0(CCIdir_img, "/", CToi, "_asS_bubble.pdf"), width = 24, height = 12)
netVisual_bubble(tPR03, targets.use = CToi)+coord_flip()
ggsave(paste0(CCIdir_img, "/", CToi, "_asT_bubble.pdf"), width = 35, height = 12)

source_path=c("THBS","GALECTIN","CXCL","LXA4","MIF")
target_path=c("SELPLG","THBS","ICAM","VWF","GALECTIN")

pdf(width = 4,height = 12,file = paste0(CCIdir_img,"/pathway_for_platelet.","as_target",".circle.05",'.pdf'))
par(mfrow = c(5,1), mar=c(0.3,0.3,1.5,0.3), xpd=TRUE)
for(ppp in target_path){
  netVisual_aggregate(tPR03, signaling = ppp,  targets.use = "Platelet", reduce = 1,
                      layout = "circle",  edge.width.max = 10, signaling.name = paste0(ppp))
}
dev.off()

### pickdotvis (only for comparison!!!)
pdf(width = 4.5,height = 9,file = paste0("./","/CM3_DotSelected_LRs_for_MphRGS1.","as_source",'.pdf'))
CellChat::netVisual_bubble(tPR03, remove.isolate = F, sources.use = c("Mph02_RGS1"),
                           targets.use = levels(tPR03@idents), #targets.use =  c("CNTN4_NEL_Cells","EC01_SELE","EC02_TLL1","Mast"),
                 pairLR.use = data.frame(interaction_name=c("Cholesterol-Cholesterol-LIPA_RORA","CXCL2_ACKR1","CXCL3_ACKR1","CXCL8_ACKR1")), angle.x = 45)+coord_flip()
CellChat::netVisual_bubble(tPR03, remove.isolate = F, sources.use = c("Mph03_FOLR2"),
                           targets.use = levels(tPR03@idents), #targets.use =  c("CNTN4_NEL_Cells","EC01_SELE","EC02_TLL1","Mast"),
                           pairLR.use = data.frame(interaction_name=c("Cholesterol-Cholesterol-LIPA_RORA","CXCL2_ACKR1","CXCL3_ACKR1","CXCL8_ACKR1")), angle.x = 45)+coord_flip()
dev.off()
pdf(width = 4.5,height = 9,file = paste0("./","/CM3_DotSelected_LRs_for_CNTN4.","as_target",'.pdf'))
CellChat::netVisual_bubble(tPR03, remove.isolate = F, #sources.use = c("B01_MS4A1","B02_IGHD","B03_MZB1","DNT_THEMIS",#"CD4T01_SELL",
                                                  ##"CNTN4_NEL_Cells" "EC01_SELE"       "EC02_TLL1"       "EC04_FBLN5"
                           sources.use = levels(tPR03@idents),
                           targets.use =  c("CNTN4_NEL_Cells"),
                           pairLR.use = data.frame(interaction_name=c("Cholesterol-Cholesterol-LIPA_RORA","CXCL2_ACKR1","CXCL3_ACKR1","CXCL8_ACKR1")), angle.x = 45)+coord_flip()
dev.off()

pdf(width = 4.5,height = 9,file = paste0("./","/CM4_DotSelected_LRs_for_CNTN4.","as_source",'.pdf'))
CellChat::netVisual_bubble(tPR03, remove.isolate = F,
                           sources.use =  c("CNTN4_NEL_Cells"),
                           targets.use = levels(tPR03@idents),
                           pairLR.use = data.frame(interaction_name=c("APP_CD74")), angle.x = 45)+coord_flip()
dev.off()
pdf(width = 6,height = 9,file = paste0("./","/CM1_DotSelected_LRs_for_CNTN4.","as_source",'.pdf'))
CellChat::netVisual_bubble(tPR03, remove.isolate = F,
                           sources.use =  c("CNTN4_NEL_Cells"),
                           targets.use = levels(tPR03@idents),
                           pairLR.use = data.frame(interaction_name=c("APP_CD74")), angle.x = 45)+coord_flip()
dev.off()
