# ============================================================================
# Disease and virus term association at cell-subtype resolution
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Association analysis of disease and virus infection related gene sets"
# Source     : extracted from Analysis240525.r (L11465-12075)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
colnames(metab_df)
newAnno_Vir=read.csv("./DiseaseAna/list2_virus.csv")
newAnno_Dis=read.csv("./DiseaseAna/list1_disease.csv")
newAnno_Vir$Location_putative="unknown"
newAnno_Vir$Ontology=newAnno_Vir$category
newAnno_Dis$Location_putative=newAnno_Dis$location
newAnno_Dis$Location_putative=gsub('""',"unknown",newAnno_Dis$Location_putativ)
newAnno_VirDis=rbind(newAnno_Vir[,c("Name","Ontology","Location_putative","System", "organ_1st","organ_2nd")],
                     newAnno_Dis[,c("Name","Ontology","Location_putative","System", "organ_1st","organ_2nd")])
#                                newAnno_VirDis$Location_putative %!in% c("TB"),]
rownames(newAnno_VirDis)=newAnno_VirDis$Name
newAnno_VirDis$Name=NULL

cor_xxx_plot=pheatmap(cor_xxx, cluster_rows = T, cluster_cols = T,
                      cutree_rows = 30, cutree_cols = 30,
         show_rownames = F, show_colnames = F,
         annotation_row = newAnno_VirDis[rownames(newAnno_VirDis) %in% colnames(cor_xxx),],
         breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
         color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100)
          )

cor_xxx_plot_clu <- row_order(cor_xxx_plot)
#  out
comb_df=data.frame()
for(i in 1:length(cor_xxx_plot_clu)){
  out <- data.frame(Items = rownames(cor_xxx[cor_xxx_plot_clu[[i]],]),
                    Cluster = rep(paste0("cluster", i), length(rownames(cor_xxx[cor_xxx_plot_clu[[i]],]))), stringsAsFactors = FALSE)
  comb_df=rbind(comb_df,out)
}
head(comb_df)

newAnno_VirDis_bak=newAnno_VirDis
newAnno_VirDis_bak$Items=rownames(newAnno_VirDis_bak)
comb_df=left_join(comb_df, newAnno_VirDis_bak, by="Items")

comb_df$Onto_clu=paste0(comb_df$Ontology, "_",comb_df$Cluster)
comb_df_bak=comb_df
rownames(comb_df_bak)=comb_df_bak$Items
comb_df_bak$Onto_clu=NULL
comb_df_bak$Items=NULL

pool_update=c()
for(j in names(table(comb_df$Onto_clu))){

  j1_onto=gsub("_.+","",j)
  j2_clu=gsub(".+_","",j)

  comb_df_tmp=comb_df[ comb_df$Ontology %in% c(j1_onto) &
                       comb_df$Cluster %in% c(j2_clu),]
  if(length(comb_df_tmp$Items)>=2){
    tmp_mtx=cor_xxx[comb_df_tmp$Items,comb_df_tmp$Items]
    tmp_mtx_vec=as.vector(tmp_mtx)
    tmp_mtx_mean=mean(tmp_mtx_vec[tmp_mtx_vec!=1])
    if(tmp_mtx_mean>=0.5){
      xxx=comb_df_tmp$Items[grepl("virus",comb_df_tmp$Items)]
      yyy=comb_df_tmp$Items[!grepl("virus",comb_df_tmp$Items)]

      pool_update=c(pool_update,xxx[1],yyy)
    }else{
      pool_update=c(pool_update,comb_df_tmp$Items)
    }
  }else{
  }

}

print(pool_update)
cor_xxx_plot_update=pheatmap(cor_xxx[pool_update,pool_update],
                      cluster_rows = T, cluster_cols = T,
                      show_rownames = T, show_colnames = F,
                      annotation_row = comb_df_bak[rownames(comb_df_bak) %in% pool_update, ],
                      breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
                      color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100))
runFig(paste0("./FigureXXX", "_clu.label"), cor_xxx_plot_update, 25, 35)

cor_xxx_plot2 <- row_order(cor_xxx_plot_update)
cor2_reorder=pool_update[cor_xxx_plot2]

target_vecs=c(cor2_reorder[1:6],
              cor2_reorder[22:29],
              cor2_reorder[43], #[40:49],
              cor2_reorder[c(57:59,63:65,69)],
              cor2_reorder[c(93,108,109,110,121:124)],
              cor2_reorder[174:180],
              cor2_reorder[208:210],
              cor2_reorder[c(248,267,275:277)])

target_vecs=c("Severe acute respiratory syndrome-related coronavirus",
   "SARS coronavirus Tor2",
   "Human coronavirus 229E",
   "Human coronavirus NL63",
   "human adenovirus 52",
   "Hepatitis C virus subtype 2k",
   "dengue virus type 1",
   "Rhinovirus A81",
   "Human bocavirus",
   "Human Respiratory syncytial virus 9320",
   "Human immunodeficiency virus 1",
   "Measles virus genotype D4",
   "Human alphaherpesvirus 3",
   "Human betaherpesvirus 7",
   "Hepatitis, Chronic",
   "Cholangitis, Sclerosing",
   "Hepatitis, Autoimmune",
   "Giant Cell Arteritis",
   "Pneumonia, Pneumocystis",
   "Adenomyosis",
   "Lupus Erythematosus, Cutaneous",
   "Lupus Nephritis",
   "Cholangiocarcinoma",
   "Pancreatic Neoplasms",
   "Lung Neoplasms",
   "Inflammatory Bowel Diseases",
   "Colorectal Neoplasms",
   "Adenomatous Polyposis Coli",
   "Urinary Bladder Neoplasms",
   "Idiopathic Pulmonary Fibrosis",
   "Stomach Neoplasms",
   "Gastroparesis",
   "Pulmonary Embolism",
   "Pneumonia, Bacterial",
   "Pancreatitis, Acute Necrotizing",
   "Gastritis, Atrophic",
   "Cystitis",
   "Skin Diseases, Infectious",
   "Enteritis",
   "Aortic Aneurysm, Abdominal",
   "Acute Coronary Syndrome",
   "Aortic Aneurysm",
   "Aortic Aneurysm, Thoracic",
   "Takayasu Arteritis",
   "Coronary Aneurysm")
target_vecs=unique(c(target_vecs, diseaseOI01$Trait))

tmp_anno=comb_df_bak[rownames(comb_df_bak) %in% target_vecs, c("System","Location_putative")] #c("Ontology","System","Location_putative")]
addAnno_vec=c("Corona virus infection",
              "Virus infection",
              "Auto-immune related disease",
              "Cancer and precancerous lesions",
              "Inflammation related disease",
              "Hemodynamic disorders")
my_Ontology=c(rep(addAnno_vec[1],4),
                    rep(addAnno_vec[2],10),
                    rep(addAnno_vec[3],8),
                    rep(addAnno_vec[4],9),
                    rep(addAnno_vec[5],8),
                    rep(addAnno_vec[6],6))
names(my_Ontology)=target_vecs[row_order(cor_xxx_plot_update2) %>% unlist()]
tmp_anno$Ontology=as.vector(my_Ontology[rownames(tmp_anno)])

cor_xxx_plot_update2=pheatmap(cor_xxx[target_vecs,target_vecs], name="rho",border_color = "white",
                             cluster_rows = T, cluster_cols = T,
                             cutree_rows = 6, cutree_cols = 6,
                             clustering_method = "ward.D2",
                             show_rownames = T, show_colnames = F,
                             annotation_row = tmp_anno,
                             breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
                             color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100))
runFig(paste0("./FigureXXX", ".cor_selected.label"), cor_xxx_plot_update2, 12, 7.5)

########### 20250527
diseaseOI01_V2=diseaseOI01
addItems=c("Atherosclerosis","Arteriosclerosis",
           "Hyperuricemia","Cholelithiasis","Cholecystolithiasis",
           "Diabetes Mellitus, Type 2","Obesity","Hypertension",
           "Hyperuricemia","Hypoglycemia","Lipid Metabolism Disorders",
           "Insulin Resistance")
addItems_anno00=data.frame(row.names=c("Atherosclerosis","Arteriosclerosis"),
                           System=c("CS","CS"),
                           Location_putative=c("CS","CS"),
                           Ontology=c("Hemodynamic disorders","Hemodynamic disorders"))
addItems_anno=diseaseOI01_V2[diseaseOI01_V2$Name %in% addItems,]
rownames(addItems_anno)=addItems_anno$Name
addItems_anno$Trait="DS"
addItems_anno$Name=NULL
colnames(addItems_anno)=colnames(addItems_anno00)

addAnno_vec=c("Corona virus infection",
              "Virus infection",
              "Auto-immune related disease",
              "Cancer and precancerous lesions",
              "Inflammation related disease",
              "Hemodynamic disorders")
my_Ontology=c(rep(addAnno_vec[1],4),
              rep(addAnno_vec[2],10),
              rep(addAnno_vec[3],8),
              rep(addAnno_vec[4],9),
              rep(addAnno_vec[5],8),
              rep(addAnno_vec[6],6))
target_vecs=c("Severe acute respiratory syndrome-related coronavirus",
              "SARS coronavirus Tor2",
              "Human coronavirus 229E",
              "Human coronavirus NL63",
              "human adenovirus 52",
              "Hepatitis C virus subtype 2k",
              "dengue virus type 1",
              "Rhinovirus A81",
              "Human bocavirus",
              "Human Respiratory syncytial virus 9320",
              "Human immunodeficiency virus 1",
              "Measles virus genotype D4",
              "Human alphaherpesvirus 3",
              "Human betaherpesvirus 7",
              "Hepatitis, Chronic",
              "Cholangitis, Sclerosing",
              "Hepatitis, Autoimmune",
              "Giant Cell Arteritis",
              "Pneumonia, Pneumocystis",
              "Adenomyosis",
              "Lupus Erythematosus, Cutaneous",
              "Lupus Nephritis",
              "Cholangiocarcinoma",
              "Pancreatic Neoplasms",
              "Lung Neoplasms",
              "Inflammatory Bowel Diseases",
              "Colorectal Neoplasms",
              "Adenomatous Polyposis Coli",
              "Urinary Bladder Neoplasms",
              "Idiopathic Pulmonary Fibrosis",
              "Stomach Neoplasms",
              "Gastroparesis",
              "Pulmonary Embolism",
              "Pneumonia, Bacterial",
              "Pancreatitis, Acute Necrotizing",
              "Gastritis, Atrophic",
              "Cystitis",
              "Skin Diseases, Infectious",
              "Enteritis",
              "Aortic Aneurysm, Abdominal",
              "Acute Coronary Syndrome",
              "Aortic Aneurysm",
              "Aortic Aneurysm, Thoracic",
              "Takayasu Arteritis",
              "Coronary Aneurysm")
names(my_Ontology)=target_vecs
tmp_annoV2=comb_df_bak[rownames(comb_df_bak) %in% target_vecs, c("System","Location_putative")] #c("Ontology","System","Location_putative")]
tmp_annoV2$Ontology=as.vector(my_Ontology[rownames(tmp_annoV2)])

tmp_annoV2=rbind(tmp_annoV2,addItems_anno00,addItems_anno)

target_vecs_V2=c(target_vecs,
                 rownames(addItems_anno00),
                 rownames(addItems_anno))

VirDis_filter_cor=read.xlsx("./VirDis_enrich_addMetab.cor0423.xlsx")
VirDis_filter_cor=VirDis_filter_cor[!is.na(VirDis_filter_cor$cor_rho),]
VirDis_filter_cor_mat=tidyr::pivot_wider(VirDis_filter_cor[!duplicated(paste0(VirDis_filter_cor$Name,VirDis_filter_cor$CT_Loc)),][,c("Name","CT_Loc","cor_rho")], names_from = Name, values_from = cor_rho)
VirDis_filter_cor_mat_final=VirDis_filter_cor_mat[,-1] %>% as.matrix()
rownames(VirDis_filter_cor_mat_final)=VirDis_filter_cor_mat$CT_Loc
VirDis_filter_cor_mat_final[is.na(VirDis_filter_cor_mat_final)]=0
cor_xxx=cor(VirDis_filter_cor_mat_final)
cor_xxx[is.na(cor_xxx)]=0

target_vecs_V3=target_vecs_V2
target_vecs_V3=c(
  target_vecs_V3[1:4],
  target_vecs_V3[5:14],
  target_vecs_V3[32:39],
  target_vecs_V3[c(22,15:21)],
  target_vecs_V3[c(40,42:43,41,44:45)],
  target_vecs_V3[c(46,48:50)],
  target_vecs_V3[23:31]
)
tmp_annoV3=tmp_annoV2
tmp_annoV3[tmp_annoV3$Ontology %in% c("metabolism"),]$Ontology="Metabolism"
tmp_annoV3[rownames(tmp_annoV3) %in% c("Atherosclerosis"),]$Ontology="Metabolism"
cor_xxx_plot_update3=pheatmap(cor_xxx[target_vecs_V3,target_vecs_V3], name="rho",border_color = "white",
                              cluster_rows = F, cluster_cols = F,
                              clustering_method = "ward.D2",
                              show_rownames = T, show_colnames = F,
                              annotation_row = tmp_annoV3[target_vecs_V3,], ###[rownames(tmp_annoV3) %in% target_vecs_V3,],
                              breaks = c(-0.5, -0.1, 0, 0.1, 0.5),
                              color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100))
runFig(paste0("./FigureXXX", ".cor_selected.label0627"), cor_xxx_plot_update3, 12, 7.5)

target_vecs=target_vecs_V3
VirDis_filter_cor_bak=VirDis_filter_cor
VirDis_filter_cor_bak$cor_rho_original=VirDis_filter_cor_bak$cor_rho
VirDis_filter_cor_bak$cor_rho=VirDis_filter_cor_bak$cor_rho*VirDis_filter_cor_bak$cor_rho
VirDis_filter_cor_bak$cor_rho=rescale(VirDis_filter_cor_bak$cor_rho, to = c(0,1))
cor_xxx3 <- list(1:4,
                 5:14,
                 15:22,
                 23:30,
                 31:36,
                 37:40,
                 41:49
                 )
comb_df3=data.frame()
for(i in 1:length(cor_xxx3)){
  out <- data.frame(Items = target_vecs[cor_xxx3[[i]]],
                    Cluster = rep(paste0("cluster", i), length(target_vecs[cor_xxx3[[i]]])), stringsAsFactors = FALSE)
  comb_df3=rbind(comb_df3,out)
}
head(comb_df3)

cor_short=VirDis_filter_cor_bak[VirDis_filter_cor_bak$cor_pval<=0.05 & VirDis_filter_cor_bak$cor_rho>=0.01,]
cor_short=cor_short[!is.na(cor_short$Name),]

#  #i="cluster1"
#  xxx=cor_short[cor_short$Name %in% alltmp_trait,]
#  xxx$GeneRatio=lapply(strsplit(xxx$Genes,split = ","), function(i){length(i)/length(allgenes_vec)}) %>% unlist()
#  xxx$cluster=i
#
#}
poolDF01=data.frame()
for(i in names(table(comb_df3$Cluster))){
  alltmp_trait=comb_df3[comb_df3$Cluster %in% i,]$Items
  xxx=cor_short[cor_short$Name %in% alltmp_trait,]
  allgenes_vec=strsplit(xxx$Genes,split = ",") %>% unlist() %>% unique()
  #xxx$GeneRatio=lapply(strsplit(xxx$Genes,split = ","), function(i){length(i)/length(allgenes_vec)}) %>% unlist()
  #xxx$cluster=i
  ct_V=unique(xxx$Celltype)
  cl_V=rep(i, length(ct_V))
  gr_V=c()
  for(x in ct_V){
    inTmp=xxx[xxx$Celltype %in% x,]
    inTmp_v=strsplit(inTmp$Genes,split = ",") %>% unlist() %>% unique()
    inTmp_Rate=length(inTmp_v)/length(allgenes_vec)
    gr_V=c(gr_V, inTmp_Rate)

    inTmp$GeneRatio=lapply(strsplit(inTmp$Genes,split = ","), function(j){length(j)/length(inTmp_v)}) %>% unlist()
    inTmp$cluster=i

  }
  com_DF=data.frame(
    cluster=cl_V,
    celltype=ct_V,
    ratio=gr_V
  )
  poolDF01=rbind(poolDF01,com_DF)
}

poolDF_tr=poolDF01[,c("celltype","ratio","cluster")] %>% group_by(celltype, cluster) %>% summarise_at("ratio",sum)

ImmuneCs=unique(aaa[aaa$MajorClass02 %in% c("Immune Cells"),]$subtypes015) %>% as.vector()
EpitheCs=unique(aaa[aaa$MajorClass02 %in% c("Epithelial Cells"),]$subtypes015) %>% as.vector()
StromaCs=unique(aaa[aaa$MajorClass02 %in% c("Stromal Cells","Endothelial Cells"),]$subtypes015) %>% as.vector()

xxxdf=poolDF_tr[poolDF_tr$celltype %in% NonImmuneCs,] %>% group_by(cluster) %>% top_n(20)
xxxdf=xxxdf %>%
  # 1. Remove grouping
  ungroup() %>%
  # 2. Arrange by
  #   i.  facet group (tissue)
  #   ii. value (score)
  arrange(cluster, -ratio) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())
ppp=ggplot(xxxdf)+
  geom_bar(aes(order, ratio), stat = "identity", fill="grey", color="white")+
  scale_x_continuous(
    breaks = xxxdf$order,
    labels = gsub("-.+","",xxxdf$celltype))+
  facet_wrap(~cluster, ncol = 4, scales = "free_x")+
  theme_pubr(border = T)+
  xlab("")+ylab("Enrichment ratio")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=10))

write.xlsx(poolDF_tr[poolDF_tr$celltype %in% ImmuneCs,],
           file = "../upload_oss/uploadRes_250627/Immune.DisVirEnrichforCelltypes.xlsx")
write.xlsx(poolDF_tr[poolDF_tr$celltype %in% NonImmuneCs,],
           file = "../upload_oss/uploadRes_250627/NonImmune.DisVirEnrichforCelltypes.xlsx")

xxxdf_addDot=xxxdf[,c("order", "cluster", "CT_Loc")]
xxxdf_addDot$group=gsub(".+-","",xxxdf_addDot$CT_Loc)
qqq=ggplot(xxxdf_addDot)+geom_point(aes(order,0, color=group), size=5)+
  scale_x_continuous(
    breaks = xxxdf_addDot$order,
    labels = gsub("-.+","",xxxdf_addDot$CT_Loc))+
  facet_wrap(~cluster, ncol = 6, scales = "free_x")+
  theme_nothing()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

mylegend=cowplot::get_legend(ggplot(xxxdf_addDot)+geom_point(aes(order,0, color=group), size=5)+
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.position = "right",
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white"),
        legend.background = element_blank(), #element_rect(fill = "white"),
        legend.title = element_blank()))

ppqq=plot_grid(ppp+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
          qqq+theme(axis.text.x = element_text(size=10)),
          align = "v", axis = "lr", ncol = 1, rel_heights = c(6,3))
ppqq_leg=plot_grid(ppqq, mylegend, ncol = 2, align = "h", axis = "none", rel_widths = c(8,1.5))

### 250821
poolDF_tr1=read.xlsx("../upload_oss/uploadRes_250627/Immune.DisVirEnrichforCelltypes.xlsx")
poolDF_tr2=read.xlsx("../upload_oss/uploadRes_250627/NonImmune.DisVirEnrichforCelltypes.xlsx")
xxxdf=poolDF_tr1[poolDF_tr1$celltype %in% ImmuneCs,] %>% group_by(cluster) %>% top_n(10)
xxxdf_Imm=xxxdf %>% ungroup() %>% arrange(cluster, -ratio) %>% mutate(order = row_number())
xxxdf=poolDF_tr2[poolDF_tr2$celltype %in% EpitheCs,] %>% group_by(cluster) %>% top_n(10)
xxxdf_Epi=xxxdf %>% ungroup() %>% arrange(cluster, -ratio) %>% mutate(order = row_number())
xxxdf=poolDF_tr2[poolDF_tr2$celltype %in% StromaCs,] %>% group_by(cluster) %>% top_n(10)
xxxdf_Str=xxxdf %>% ungroup() %>% arrange(cluster, -ratio) %>% mutate(order = row_number())

poolDF01=data.frame()
poolDF02=data.frame()
for(i in names(table(comb_df3$Cluster))){
  alltmp_trait=comb_df3[comb_df3$Cluster %in% i,]$Items
  xxx=cor_short[cor_short$Name %in% alltmp_trait,]
  xxx=xxx[xxx$Celltype %in% ImmuneCs,]
  #xxx=xxx[xxx$Celltype %in% EpitheCs,]
  #xxx=xxx[xxx$Celltype %in% StromaCs,]
  allgenes_vec=strsplit(xxx$Genes,split = ",") %>% unlist() %>% unique()

  ct_V=unique(xxx$Celltype)
  cl_V=rep(i, length(ct_V))
  gr_V=c()
  for(x in ct_V){
    inTmp=xxx[xxx$Celltype %in% x,]
    inTmp_v=strsplit(inTmp$Genes,split = ",") %>% unlist() %>% unique()
    inTmp_Rate=length(inTmp_v)/length(allgenes_vec)
    gr_V=c(gr_V, inTmp_Rate)

    inTmp$GeneRatio=lapply(strsplit(inTmp$Genes,split = ","), function(j){length(j)/length(inTmp_v)}) %>% unlist()
    inTmp$cluster=i

  }
  com_DF=data.frame(
    cluster=cl_V,
    celltype=ct_V,
    ratio=gr_V
  )
  poolDF01=rbind(poolDF01,com_DF)

  loc_V=unique(xxx$Location)
  cl_V0=rep(i, length(loc_V))
  gr_V0=c()
  for(x in loc_V){
    inTmp=xxx[xxx$Location %in% x,]
    inTmp_v=strsplit(inTmp$Genes,split = ",") %>% unlist() %>% unique()
    inTmp_Rate=length(inTmp_v)/length(allgenes_vec)
    gr_V0=c(gr_V0, inTmp_Rate)

    inTmp$GeneRatio=lapply(strsplit(inTmp$Genes,split = ","), function(j){length(j)/length(inTmp_v)}) %>% unlist()
    inTmp$cluster=i

  }
  com_DF=data.frame(
    cluster=cl_V0,
    location=loc_V,
    ratio=gr_V0
  )
  poolDF02=rbind(poolDF02,com_DF)

}

poolDF_tr=poolDF[,c("CT_Loc","GeneRatio","cluster")] %>% group_by(CT_Loc, cluster) %>% summarise_at("GeneRatio",sum)

celltypes_list=list(ImmuneCs=ImmuneCs,
                    EpitheCs=EpitheCs,
                    StromaCs=StromaCs)

xxxdf=poolDF_tr[gsub("-.+","",poolDF_tr$CT_Loc) %in% StromaCs,] %>% group_by(cluster) %>% top_n(20)
xxxdf=xxxdf %>%
  # 1. Remove grouping
  ungroup() %>%
  # 2. Arrange by
  #   i.  facet group (tissue)
  #   ii. value (score)
  arrange(cluster, -ratio) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())
ppp=ggplot(xxxdf)+
  geom_bar(aes(order, ratio), stat = "identity", fill="grey", color="white")+
  scale_x_continuous(
    breaks = xxxdf$order,
    labels = gsub("-.+","",xxxdf$celltype))+
  facet_wrap(~cluster, ncol = 4, scales = "free_x")+
  theme_pubr(border = T)+
  xlab("")+ylab("Enrichment ratio")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=10))

xxxdf_addDot=xxxdf[,c("order", "cluster", "CT_Loc")]
xxxdf_addDot$group=gsub(".+-","",xxxdf_addDot$CT_Loc)
qqq=ggplot(xxxdf_addDot)+geom_point(aes(order,0, color=group), size=5)+
  scale_x_continuous(
    breaks = xxxdf_addDot$order,
    labels = gsub("-.+","",xxxdf_addDot$CT_Loc))+
  facet_wrap(~cluster, ncol = 6, scales = "free_x")+
  theme_nothing()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

mylegend=cowplot::get_legend(ggplot(xxxdf_addDot)+geom_point(aes(order,0, color=group), size=5)+
                               guides(colour = guide_legend(override.aes = list(shape = 15))) +
                               theme(legend.position = "right",
                                     legend.key = element_blank(),
                                     legend.box.background = element_rect(fill = "white"),
                                     legend.background = element_blank(), #element_rect(fill = "white"),
                                     legend.title = element_blank()))

ppqq=plot_grid(ppp+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
               qqq+theme(axis.text.x = element_text(size=10)),
               align = "v", axis = "lr", ncol = 1, rel_heights = c(6,3))
ppqq_leg=plot_grid(ppqq, mylegend, ncol = 2, align = "h", axis = "none", rel_widths = c(8,1.5))

data_long=xxxdf_Imm
colnames(data_long)=c("source","target","value")
data_long$value=rescale(data_long$value-median(data_long$value), to = c(0.01,1.01))
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
data_long$IDsource=match(data_long$source, nodes$name)-1
data_long$IDtarget=match(data_long$target, nodes$name)-1
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
p=sankeyNetwork(Links = data_long, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=20, nodePadding=20)
## Make a webshot in pdf : high quality but can not choose printed zone
saveNetwork(p,"./tmp00.html")
webshot("./tmp00.html", "./CTweight_ImmC.pdf") #, cliprect = c(0, 100, 1000, 10))

#  1 / (1 + exp(-x))
#}
