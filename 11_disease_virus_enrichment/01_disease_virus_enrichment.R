# ============================================================================
# DisGeNET and virus-infection risk gene set overlap with hypergeometric enrichment
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Association analysis of disease and virus infection related gene sets"
# Source     : extracted from Analysis240525.r (L9204-9339, L11147-11463)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
setwd("./DiseaseAna/")

AllDis=read.table("./Disgenenet_filtered.txt", sep="\t", header = T, fill = TRUE)

VisrusDF=read.csv("./phipster_LR100.csv")
VisrusDF$VirusID=gsub("\\|.+","",VisrusDF$interaction_pair)
VisrusDF$Gene=gsub(".+\\|","",VisrusDF$interaction_pair)
VisrusDF_tmp=VisrusDF[,c("VirusID","Gene")]
VisrusDF_tmp$Gene=gsub(".+ ","",VisrusDF_tmp$Gene)
VisrusDF_tmp=unique(VisrusDF_tmp)

#disease enrich
D2gene=VisrusDF_tmp
D2gene$ID=D2gene$VirusID

D2gene=AllDis
D2gene=D2gene[,c("geneSymbol","diseaseName")] %>% unique()
D2gene$ID=D2gene$diseaseName
D2gene$Gene=D2gene$geneSymbol

DEGlistFile_tmp=DEGs_Mobj_33CT
DEGlistFile=DEGlistFile_tmp
total_background=34016

dotdf=data.frame()
Trait=c()
Celltype=c()
Overlap=c()
Pvalue=c()
OverGenes=c()
for (d in 1:length(unique(D2gene$ID))){
  #d=1
  myDname=unique(D2gene$ID)[d]
  tmpdf_snp=unique(D2gene[D2gene$ID==myDname,]$Gene)
  tmpdf_snp=tmpdf_snp[!is.na(tmpdf_snp)]

  #  #c="FB_APOC1"
  for(c in unique(DEGlistFile_tmp$cluster)){
    myDEGlist=DEGlistFile_tmp[DEGlistFile_tmp$cluster==c,]
        df=myDEGlist
        if (nrow(df) < 1) next
        df <- df[df$p_val_adj < 0.1, ]
        myDEGlist <-  df
        if (nrow(myDEGlist) < 1) next

        goverlap=intersect(tmpdf_snp, myDEGlist$gene)
        pvalue=phyper(length(goverlap)-1, length(myDEGlist$gene)-length(goverlap),
                      total_background-length(myDEGlist$gene), length(tmpdf_snp),lower.tail = F)

        Trait=c(Trait, myDname)
        Celltype=c(Celltype, c)
        Overlap=c(Overlap, length(goverlap))
        Pvalue=c(Pvalue, pvalue)
        ggg=if_else(length(goverlap)>0, paste(goverlap, sep="", collapse = ","), "NA")
        OverGenes=c(OverGenes, ggg)
    #  }
    #}
  }
}

dotdf=data.frame(Trait=Trait,  Celltype=Celltype,
                 Overlap=Overlap, Pvalue=Pvalue, Genes=OverGenes)
dotdf$`-log10(FDR)`=-log10(qvalue::qvalue(p = dotdf$Pvalue, pi0 = 1)$qvalues)
dotdf$`-log10(Pvalue)`=-log10(dotdf$Pvalue)
dotdf$`#DEGs\nenriched`=dotdf$Overlap

dotdf=readRDS("./DiseaseAna/Disgenenet_filtered.enrich.rds")
xxx=dotdf[dotdf$`-log10(FDR)`>=-log10(0.05),]

ggplot(xxx[xxx$Trait %in% c("Adenocarcinoma", "Breast Neoplasms", "Hepatolenticular Degeneration"),])+geom_point(aes(Celltype, Trait, fill=`#DEGs\nenriched`, size=`-log10(Pvalue)`), shape=21)+
  theme_minimal()+
  scale_size_continuous(limits = c(0.01,NA))+
  scale_fill_gradient(low = "navy", high = "yellow")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size=14))+
  theme(axis.text.y = element_text(colour = "black", size=14))+
  theme(strip.text = element_text(colour = "black", size=14))+
  xlab("")+ylab("")#+

xxx=xxx[xxx$Overlap>3,]
pval_vec=c()
rho_vec=c()
for(x in 1:dim(xxx)[1]){
  mygenes=strsplit(xxx[x,]$Genes,",") %>% unlist()
  tDF=AllDis[AllDis$diseaseName %in% c(xxx[x,]$Trait) & AllDis$geneSymbol %in% mygenes,]
  cDF=DEGlistFile_tmp[DEGlistFile_tmp$cluster %in% c(xxx[x,]$Celltype) & DEGlistFile_tmp$gene %in% mygenes,]
  tcDF=left_join(tDF,cDF, by=c("geneSymbol"="gene"))
  mycorr=cor.test(tcDF$score, tcDF$avg_log2FC, method = "spearman")
  pval_vec=c(pval_vec, as.vector(mycorr$p.value))
  rho_vec=c(rho_vec, as.vector(mycorr$estimate))
}
xxx$cor_pval=pval_vec
xxx$cor_rho=rho_vec
xxx_filter=xxx[xxx$cor_pval<=0.05 & xxx$cor_rho>=0.25,]
xxx_filter=xxx_filter[!is.na(xxx_filter$Trait),]
xxx_filter=xxx_filter[order(xxx_filter$Celltype, xxx_filter$cor_rho, decreasing = c(FALSE, TRUE), method=c("radix")), ]
xxx_filter$Trait=factor(xxx_filter$Trait, levels = unique(xxx_filter$Trait))
ggplot(xxx_filter)+geom_point(aes(cor_rho, Trait, size=cor_rho, color=-log(cor_pval), fill=-log(cor_pval)))+
  facet_wrap(~Celltype, ncol = 8, scales = "free")
ggsave(filename = "./Disgenenet_filtered.enrich.cor_add.pdf",width = 40, height = 24)

disgen_filter=readRDS("./DiseaseAna/Disgenenet_filtered.enrich.rds")
phipster_filter=readRDS("./DiseaseAna/phipster_LR100.enrich.rds")

VisrusDF=read.csv("./DiseaseAna/phipster_LR100.csv")
VisrusDF$VirusID=gsub("\\|.+","",VisrusDF$interaction_pair)
VisrusDF$Gene=gsub(".+\\|","",VisrusDF$interaction_pair)
VisrusDF$Gene=gsub(".+ ","",VisrusDF$Gene)
VisrusDF=unique(VisrusDF)
VisrusDF=VisrusDF[,c("VirusID","Gene","Final_LR")]
VisrusDF$Type="Virus"
colnames(VisrusDF)=c("trait","gene","score","type")

### many items in vnames have NO common names; actually taxid refers to NCBI taxonomy id!!! 20250313

VisrusDF$oriID=VisrusDF$trait
VisrusDF$trait=gsub("_.+","",VisrusDF$trait)

AllDisDF=AllDis[,c("diseaseName","geneSymbol","score","X","diseaseId")]
colnames(AllDisDF)=c("trait","gene","score","type","oriID")
AllDisDF$type="Disease"

D2gene=rbind(VisrusDF, AllDisDF)
D2gene$ID=D2gene$trait
D2gene$Gene=D2gene$gene

DEGlistFile_tmp=DEGs_Mobj_107subCT #DEGs_Mobj_33CT

total_background=34016  ### unique(DEGs_Mobj_107subCT$gene) %>% length():  16923

dotdf=data.frame()
Trait=c()
Celltype=c()
Location=c()
Overlap=c()
Pvalue=c()
OverGenes=c()
for (d in 1:length(unique(D2gene$ID))){
  #d=1
  myDname=unique(D2gene$ID)[d]
  tmpdf_snp=unique(D2gene[D2gene$ID==myDname,]$Gene)
  tmpdf_snp=tmpdf_snp[!is.na(tmpdf_snp)]
  if (length(tmpdf_snp) < 2) next

  for(gn in unique(DEGlistFile_tmp$location)){
      myDEGlist_loc= DEGlistFile_tmp[DEGlistFile_tmp$location==gn,]
      if (nrow(myDEGlist_loc) < 3) next

  for(c in unique(myDEGlist_loc$subtypes)){
    myDEGlist=myDEGlist_loc[myDEGlist_loc$subtypes==c,]
    df=myDEGlist
    if (nrow(df) < 1) next

    df <- df[df$p_val_adj < 0.1, ]
    myDEGlist <-  df
    if (nrow(myDEGlist) < 1) next

    goverlap=intersect(tmpdf_snp, myDEGlist$gene)
    if (length(goverlap) <1) next

    pvalue=phyper(length(goverlap)-1, length(myDEGlist$gene)-length(goverlap),
                  total_background-length(myDEGlist$gene), length(tmpdf_snp),lower.tail = F)

    Trait=c(Trait, myDname)
    Celltype=c(Celltype, c)
    Location=c(Location, gn)
    Overlap=c(Overlap, length(goverlap))
    Pvalue=c(Pvalue, pvalue)
    ggg=if_else(length(goverlap)>0, paste(goverlap, sep="", collapse = ","), "NA")
    OverGenes=c(OverGenes, ggg)

    }
  }
}

dotdf=data.frame(Trait=Trait,  Celltype=Celltype,
                 Location=Location,
                 Overlap=Overlap, Pvalue=Pvalue, Genes=OverGenes)
dotdf$`-log10(FDR)`=-log10(qvalue::qvalue(p = dotdf$Pvalue, pi0 = 1)$qvalues)
dotdf$`-log10(Pvalue)`=-log10(dotdf$Pvalue)
dotdf$`#DEGs\nenriched`=dotdf$Overlap

saveRDS(dotdf,"./VirDisGenenet_raw.enrich.rds")

#xxx=dotdf[dotdf$`-log10(FDR)`>=-log10(0.05),]

ggplot(xxx[xxx$Trait %in% c("Adenocarcinoma", "Breast Neoplasms", "Hepatolenticular Degeneration"),])+geom_point(aes(Celltype, Trait, fill=`#DEGs\nenriched`, size=`-log10(Pvalue)`), shape=21)+
  theme_minimal()+
  scale_size_continuous(limits = c(0.01,NA))+
  scale_fill_gradient(low = "navy", high = "yellow")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size=14))+
  theme(axis.text.y = element_text(colour = "black", size=14))+
  theme(strip.text = element_text(colour = "black", size=14))+
  xlab("")+ylab("")#+

xxx=xxx[xxx$Overlap>3,]
pval_vec=c()
rho_vec=c()
for(x in 1:dim(xxx)[1]){
  mygenes=strsplit(xxx[x,]$Genes,",") %>% unlist()
  tDF=AllDis[AllDis$diseaseName %in% c(xxx[x,]$Trait) & AllDis$geneSymbol %in% mygenes,]
  cDF=DEGlistFile_tmp[DEGlistFile_tmp$cluster %in% c(xxx[x,]$Celltype) & DEGlistFile_tmp$gene %in% mygenes,]
  tcDF=left_join(tDF,cDF, by=c("geneSymbol"="gene"))
  mycorr=cor.test(tcDF$score, tcDF$avg_log2FC, method = "spearman")
  pval_vec=c(pval_vec, as.vector(mycorr$p.value))
  rho_vec=c(rho_vec, as.vector(mycorr$estimate))
}
xxx$cor_pval=pval_vec
xxx$cor_rho=rho_vec
xxx_filter=xxx[xxx$cor_pval<=0.05 & xxx$cor_rho>=0.25,]
xxx_filter=xxx_filter[!is.na(xxx_filter$Trait),]
xxx_filter=xxx_filter[order(xxx_filter$Celltype, xxx_filter$cor_rho, decreasing = c(FALSE, TRUE), method=c("radix")), ]
xxx_filter$Trait=factor(xxx_filter$Trait, levels = unique(xxx_filter$Trait))
ggplot(xxx_filter)+geom_point(aes(cor_rho, Trait, size=cor_rho, color=-log(cor_pval), fill=-log(cor_pval)))+
  facet_wrap(~Celltype, ncol = 8, scales = "free")
ggsave(filename = "./Disgenenet_filtered.enrich.cor_add.pdf",width = 40, height = 24)

diseaseOI=read.xlsx("./DiseaseAna/list_disease.xlsx")
diseaseOI$Trait=diseaseOI$Name
diseaseOI$Location_putative=diseaseOI$location
diseaseOI=diseaseOI[,c("Trait","Location_putative","Name","Ontology")]

diseaseOI01=metab_df
diseaseOI01$Trait=diseaseOI01$Disease
diseaseOI01$Name=diseaseOI01$Disease
diseaseOI01$Location_putative=diseaseOI01$location
diseaseOI01=diseaseOI01[,c("Trait","Location_putative","Name","Ontology")]

virusOI=read.xlsx("./DiseaseAna/list_virus.xlsx")
virusOI$Trait=paste0("tx",virusOI$ID)
virusOI$Location_putative=virusOI$potential.organ
virusOI$Ontology="virus infection"
virusOI=virusOI[,c("Trait","Location_putative","Name","Ontology")]

###20250423
VirDisGenenet_filter=readRDS("./VirDisGenenet_raw.enrich.rds")
VirDisGenenet_filter=VirDisGenenet_filter[!is.na(VirDisGenenet_filter$Pvalue),]
VirDisGenenet_filter=VirDisGenenet_filter[VirDisGenenet_filter$`-log10(FDR)`>=-log10(0.05),]
VirDisGenenet_filter=VirDisGenenet_filter[VirDisGenenet_filter$Trait %in% diseaseOI01$Trait,]
VirDisGenenet_filter[,c("Location_putative","Name","Ontology")]=NULL
diseaseOI_filter01=left_join(diseaseOI01,VirDisGenenet_filter,by = "Trait")
diseaseOI_filter01=diseaseOI_filter01[!is.na(diseaseOI_filter01$Pvalue),]

##########################################################################
#xxx=diseaseOI_filter
#xxx=virusOI_filter
xxx=diseaseOI_filter01

AllDis=D2gene

xxx=xxx[xxx$Overlap>=3,]
pval_vec=c()
rho_vec=c()
for(x in 1:dim(xxx)[1]){
  #x=3
  if(mod(x,5000)==0){
    print(x)
  }
  mygenes=strsplit(xxx[x,]$Genes,",") %>% unlist()
  if(length(mygenes)>=3){
  tDF=AllDis[AllDis$trait %in% c(xxx[x,]$Trait) & AllDis$Gene %in% mygenes,]
  cDF=DEGlistFile_tmp[DEGlistFile_tmp$subtypes %in% c(xxx[x,]$Celltype) &
                      DEGlistFile_tmp$location %in% c(xxx[x,]$Location) &
                      DEGlistFile_tmp$gene %in% mygenes,]
  tcDF=left_join(tDF,cDF, by="gene")
  mycorr=cor.test(tcDF$score, tcDF$avg_log2FC, method = "spearman")
  pval_vec=c(pval_vec, as.vector(mycorr$p.value))
  rho_vec=c(rho_vec, as.vector(mycorr$estimate))
  }
}
xxx=xxx[!is.na(xxx$Overlap),]
xxx$cor_pval=pval_vec
xxx$cor_rho=rho_vec

VirDis_filter_cor=rbind(diseaseOI_filter_cor,virusOI_filter_cor)

###20250423
VirDis_filter_cor=readRDS("./VirDis_filtered.enrich.cor_add.rds")
VirDis_filter_cor=rbind(VirDis_filter_cor, xxx)
VirDis_filter_cor$CT_Loc=paste0(VirDis_filter_cor$Celltype,"-",VirDis_filter_cor$Location)

VirDis_filter_cor=VirDis_filter_cor_bak

###20250515
VirDis_filter_cor=read.xlsx("./VirDis_enrich_addMetab.cor0423.xlsx")
VirDis_filter_cor=VirDis_filter_cor[!is.na(VirDis_filter_cor$cor_rho),]
VirDis_filter_cor=VirDis_filter_cor[VirDis_filter_cor$cor_rho>=0.2, ]

VirDis_filter_cor_mat=tidyr::pivot_wider(VirDis_filter_cor[!duplicated(paste0(VirDis_filter_cor$Name,VirDis_filter_cor$CT_Loc)),][,c("Name","CT_Loc","cor_rho")], names_from = Name, values_from = cor_rho)
VirDis_filter_cor_mat_final=VirDis_filter_cor_mat[,-1] %>% as.matrix()
rownames(VirDis_filter_cor_mat_final)=VirDis_filter_cor_mat$CT_Loc
VirDis_filter_cor_mat_final[is.na(VirDis_filter_cor_mat_final)]=0

anno_df=VirDis_filter_cor[VirDis_filter_cor$Name %in% colnames(VirDis_filter_cor_mat_final),
                          c("Name","Ontology","Location_putative")] %>% unique()
anno_df=anno_df[!duplicated(anno_df$Name),]
anno_df_in=anno_df[,-1]
rownames(anno_df_in)=anno_df$Name
anno_df_in$metab_high=anno_df_in$Ontology
anno_df_in[anno_df_in$Ontology %!in% c("metabolism"),]$metab_high="other"

cor_xxx=cor(VirDis_filter_cor_mat_final)
cor_xxx[is.na(cor_xxx)]=0
cor_xxx_plot=pheatmap(cor_xxx, cluster_rows = T, cluster_cols = T,
         show_rownames = F, show_colnames = F,
         annotation_row = anno_df_in,
         breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
         color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100)
          )
cor_xxx_reorder=cor_xxx[row_order(cor_xxx_plot),column_order(cor_xxx_plot)]
write.csv(cor_xxx_reorder, "./VirDis_cor_reorder.csv")
runFig("./xxx000", cor_xxx_plot, w = 18, h = 18)

anno_df_in$show=0
anno_df_in[rownames(anno_df_in) %in% target_vecs,]$show=1

cor_xxx_plot=pheatmap(cor_xxx, cluster_rows = T, cluster_cols = T,
                      show_rownames = T, show_colnames = F,
                      annotation_row = anno_df_in,
                      breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
                      color = colorRampPalette(c(rep("navy", 2), rep("lightblue", 1), rep("white", 2), rep("red", 2), rep("#A50000", 1)))(100)
)
runFig("./xxx000_highlight", cor_xxx_plot, w = 24, h = 96)
