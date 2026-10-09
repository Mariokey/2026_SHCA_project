# ============================================================================
# GO / KEGG enrichment of DEGs (clusterProfiler, RunHsGO, gseFun, toppEnrich)
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Function annotation"
# Source     : extracted from Analysis240525.r (L1265-1292, L1544-1591, L1914-2116, L12159-12464)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
gseFun=function(gene_list,...){

## omit any NA values
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = "org.Hs.eg.db",
             pAdjustMethod = "none")

gseDf=gse@result
gseDf=gseDf[gseDf$pvalue<=0.05,]

return(
  list(
    gse_obj=gse,
    gse_fdf=gseDf
  )
)
}

RunHsGO <- function(entrez_id,updown="up",i=0, trans=TRUE){

library(org.Hs.eg.db)
OrgDb <- org.Hs.eg.db

if(trans){
  entrez_id <- mapIds(x=OrgDb, keys = entrez_id, keytype = "SYMBOL", column = "ENTREZID")
  entrez_id <- na.omit(entrez_id)
}

if(length(entrez_id)>0){
  goresult <- enrichGO(gene = entrez_id,
                       OrgDb = OrgDb,
                       ont = "ALL",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       readable = TRUE)

  if(!is.null(goresult)){

    plotdata <- goresult@result %>% group_by(ONTOLOGY) %>% top_n(n=100, wt=-p.adjust)
    plotdata$GO_Term <- paste(plotdata$ID, plotdata$Description,sep = ' ')
    plotdata$Descript <- lapply(plotdata$GO_Term,function(x){
      x <- paste(substring(x,1,40),'...',sep='')
      x = gsub(pattern = "GO:\\d+ ",replacement = "",x,perl = TRUE)
    })
    plotdata$Descript=as.character(plotdata$Descript)

    plotdata$GO_Term <- as.character(plotdata$GO_Term)
    plotdata=plotdata[plotdata$ONTOLOGY=="BP",]
    ddd=as.data.frame(plotdata)
    ddd$Factor=sapply(ddd$GeneRatio, function(x) eval(parse(text=x)))
    pp=ggplot(ddd)+geom_point(aes(Factor, Descript, size=Count, color=pvalue))+
      scale_color_gradient(low="red", high="blue", trans="log10")+
      xlab("Rich Factor")+ylab("GO Term")
    return(list(d=ddd,p=pp))

  }
}
}

toppEnrich=function(gene_vector, maxtry=20, interSecond=5, geneflatten=FALSE){

library(httr)
library(curl)
library(jsonlite)

extract_entrez_ids <- function(x) {
  ids <- c()
  if (is.list(x)) {
    for(y in 1:length(x)){
      #y=1
      xx=x[[y]]
      if ("Entrez" %in% names(xx)) {
        ids <- c(ids, xx$Entrez)
      }
      #}
    }
  }
  return(ids)
}

tmprun1=function(json_data){output1=POST("https://toppgene.cchmc.org/API/lookup",
                                         body = json_data,
                                         add_headers('Content-Type' = 'text/json'))
return(output1)
}
#}
myrun1=function(json_data, i=1){
  if(i<=maxtry){
    tryCatch({return(tmprun1(json_data))}, error=function(e, a=i){message(paste0("try-", a));
      a=a+1;
      Sys.sleep(interSecond);
      return(myrun1(json_data, a))})
  }else{
    stop("max try!")
  }
}

tmprun2=function(json_data){output2=POST("https://toppgene.cchmc.org/API/enrich",
                                         body = json_data,
                                         add_headers('Content-Type' = 'text/json'),
                                         timeout(120))
return(output2)
}
myrun2=function(json_data, i=1){
  if(i<=maxtry){
    tryCatch({return(tmprun2(json_data))}, error=function(e, a=i){message(paste0("try-", a));
      a=a+1;
      Sys.sleep(interSecond);
      return(myrun2(json_data, a))})
  }else{
    stop("max try!")
  }
}

# Convert the gene vector to JSON format
json_data <- toJSON(list(Symbols = gene_vector))

response = myrun1(json_data)

# Check the status of the response
print(status_code(response))

# Parse and print the response content
content <- content(response, "parsed")

##############################################################################
mygenes=extract_entrez_ids(content$Genes)
print(mygenes)

annotations_flat=data.frame()

cmd01=' -H \'Content-Type: text/json\' -d \'{"Genes":['
#9956,1593,10457,8613,7941,245972,203100,284129
cmd03='], "Categories": [{"Type": "ToppCell","PValue": 0.05,"MinGenes": 1,"MaxGenes": 1500,"MaxResults": 20,"Correction": "FDR"}]}\' https://toppgene.cchmc.org/API/enrich'

outpost=system2(command = "curl", args = paste0(cmd01, paste(mygenes, collapse = ","), cmd03), stdout=TRUE)
tmpdata <- jsonlite::fromJSON(outpost, flatten = TRUE)
annotations <- tmpdata$Annotations

if(length(annotations)>0){
  annotations_flat <- annotations %>% as.data.frame()

  if(geneflatten){
    genes_expanded <- annotations_flat %>%
      mutate(row_id = row_number()) %>%
      tidyr::unnest(Genes) %>%
      select(row_id, everything()) %>% as.data.frame()
    annotations_flat=genes_expanded
  }
}

return(annotations_flat)

}

poolDF_new=data.frame()
for(i in names(table(comb_df3$Cluster))){
  alltmp_trait=comb_df3[comb_df3$Cluster %in% i,]$Items
  xxx=cor_short[cor_short$Name %in% alltmp_trait,]

###xxx=xxx[xxx$Celltype %in% NonImmuneCs,]
#xxx=xxx[xxx$Celltype %in% ImmuneCs,]
#xxx=xxx[xxx$Celltype %in% EpitheCs,]
xxx=xxx[xxx$Celltype %in% StromaCs,]

  allgenes_vec=strsplit(xxx$Genes,split = ",") %>% unlist() %>% unique()
  genes=c()
  weights=c()
  for(x in allgenes_vec){
    genes=c(genes,x)
    weights=c(weights,mean(cor_short[grepl(x,cor_short$Genes),]$cor_rho))
  }
  newDF=data.frame(gene=genes, avg_log2fc=weights)
  newDF$cluster=i
  poolDF_new=rbind(poolDF_new, newDF)
}

#step-1

#step-2
poolDF_new_ImmC_uniq=poolDF_new_ImmC
poolDF_new_ImmC_uniq$avg_log2FC=poolDF_new_ImmC_uniq$avg_log2fc
poolDF_new_ImmC_uniq = DEG_clusterUniq(poolDF_new_ImmC_uniq, minGnum = 5)
poolDF_new_EpiC_uniq=poolDF_new_EpiC
poolDF_new_EpiC_uniq$avg_log2FC=poolDF_new_EpiC_uniq$avg_log2fc
poolDF_new_EpiC_uniq = DEG_clusterUniq(poolDF_new_EpiC_uniq, minGnum = 5)
poolDF_new_StromC_uniq=poolDF_new_StromC
poolDF_new_StromC_uniq$avg_log2FC=poolDF_new_StromC_uniq$avg_log2fc
poolDF_new_StromC_uniq = DEG_clusterUniq(poolDF_new_StromC_uniq, minGnum = 5)

#step-4
poolDF_new=poolDF_new_StromC_uniq

poolDF_new_top100=poolDF_new #%>% group_by(cluster) %>% top_n(30,avg_log2fc)
poolDF_new_top100$avg_log2FC=NULL
poolDF_new_top100_funcs=enrichFUNgo(ud_table = poolDF_new_top100,intermode = FALSE,enriches = c("GO","KEGG"))

poolDF_new_StromC_top100_funcs=poolDF_new_top100_funcs
saveRDS(list(poolDF_new_ImmC_top100_funcs,
           poolDF_new_EpiC_top100_funcs,
           poolDF_new_StromC_top100_funcs), file = "./poolDF_new_ImmC_EpiC_StromC_enrichFungo.rds")

write.xlsx(list(poolDF_new_ImmC_top100_funcs$dflist$GO,
              poolDF_new_ImmC_top100_funcs$dflist$KEGG,
              poolDF_new_EpiC_top100_funcs$dflist$GO,
              poolDF_new_EpiC_top100_funcs$dflist$KEGG,
              poolDF_new_StromC_top100_funcs$dflist$GO,
              poolDF_new_StromC_top100_funcs$dflist$KEGG), file = "./ImmC_EpiC_StromC_VirDisCluster_FunEnrich.xlsx")

DEgene_list=list(
  poolDF_new_ImmC_uniq=poolDF_new_ImmC_uniq,
  poolDF_new_StromC_uniq=poolDF_new_StromC_uniq,
  poolDF_new_EpiC_uniq=poolDF_new_EpiC_uniq
)

saveRDS(list(comb_df3=comb_df3,
             cor_short=cor_short,
             celltypes_list=celltypes_list,
             DEgene_list=DEgene_list), file="./CMall4_topGene_topCell_topLoc.rds")

combined_keggDF=data.frame()

GO_funcs_sub=poolDF_new_StromC_top100_funcs$dflist$KEGG

D2gene=D2gene_bak

GO_funcs_sub$Factor=sapply(GO_funcs_sub$GeneRatio, function(x) eval(parse(text=x)))
fDF=table(GO_funcs_sub$Description,GO_funcs_sub$celltype)
fDF_names=rownames(fDF)[rowSums(fDF)<=3]
GO_funcs_sub=GO_funcs_sub[GO_funcs_sub$Description %in% fDF_names,] %>% dplyr::group_by(celltype) %>% dplyr::top_n(8,Factor) %>% as.data.frame()

unique_GO=unique(GO_funcs_sub$GO_Term)
virusOI_Trait2Name=dplyr::left_join(virusOI[,c("Trait","Name")],
                                    D2gene[,c("trait","gene")],by = c("Trait"="trait"))
virusOI_Trait2Name=virusOI_Trait2Name[,c("Name","gene")]
colnames(virusOI_Trait2Name)=c("ID","gene")
D2gene=rbind(D2gene[D2gene$type %in% c("Disease"),c("ID","gene"),],
             virusOI_Trait2Name) %>% unique()
D2gene=D2gene[D2gene$ID %in% target_vecs,]

D2gene=D2gene[D2gene$gene %in% unique(poolDF_new_ImmC_uniq$gene),]

total_background=7562

GO_Dotplot_df=data.frame()
Trait=c()
Func=c()
Overlap=c()
Enrichment=c()
Pvalue=c()
OverGenes=c()
for(mygo in unique_GO){

### GO

### KEGG
mygo_genes=strsplit(GO_funcs_sub[GO_funcs_sub$GO_Term %in% mygo,]$GeneSymbol,split = "/") %>% unlist() %>% unique()

for (d in 1:length(unique(D2gene$ID))){
  #d=1
  myDname=unique(D2gene$ID)[d]
  tmpdf_snp=unique(D2gene[D2gene$ID==myDname,]$gene)
  tmpdf_snp=tmpdf_snp[!is.na(tmpdf_snp)]
  if (length(tmpdf_snp) < 2) next

  goverlap=intersect(tmpdf_snp, mygo_genes)
  if (length(goverlap) <1) next

  pvalue=phyper(length(goverlap)-1, length(mygo_genes)-length(goverlap),
                total_background-length(mygo_genes), length(tmpdf_snp),lower.tail = F)

  Trait=c(Trait, myDname)
  Func=c(Func, mygo)
  Overlap=c(Overlap, length(goverlap))
  Enrichment=c(Enrichment, round(length(goverlap)/length(unique(c(mygo_genes,tmpdf_snp))),digits = 4)*100)
  Pvalue=c(Pvalue, pvalue)
  ggg=dplyr::if_else(length(goverlap)>0, paste(goverlap, sep="", collapse = ","), "NA")
  OverGenes=c(OverGenes, ggg)

}
}
GO_Dotplot_df=data.frame(Trait=Trait,
                       Func=Func,
                       Overlap=Overlap, Enrichment=Enrichment, Pvalue=Pvalue, Genes=OverGenes)
GO_Dotplot_df$`-log10(FDR)`=-log10(qvalue::qvalue(p = GO_Dotplot_df$Pvalue, pi0 = 1)$qvalues)
GO_Dotplot_df$`-log10(Pvalue)`=-log10(GO_Dotplot_df$Pvalue)
GO_Dotplot_df$`#Genes`=GO_Dotplot_df$Overlap
GO_Dotplot_df$Trait=factor(GO_Dotplot_df$Trait, levels = rev(target_vecs))
GO_Dotplot_df$Func=factor(GO_Dotplot_df$Func, levels = unique_GO)
GO_Dotplot_df$Enrichment=rescale(GO_Dotplot_df$Enrichment, to = c(0,1))

combined_keggDF=rbind(combined_keggDF, GO_Dotplot_df)

#GO_Dotplot_df[is.infinite(GO_Dotplot_df$Pvalue),]$Pvalue
combined_keggDF=combined_keggDF %>% group_by(Func, Parts) %>% top_n(10, Enrichment)
combined_keggDF$Parts=factor(combined_keggDF$Parts, levels=c("Immune Cells","Epithelial Cells","Stroma Cells"))

combined_keggDF_tmp=combined_keggDF
combined_keggDF_tmp$Trait=factor(combined_keggDF_tmp$Trait, levels = target_vecs)
combined_keggDF_tmp=combined_keggDF_tmp[order(combined_keggDF_tmp$Trait),]

combined_keggDF$Func=factor(combined_keggDF$Func, levels = unique(combined_keggDF_tmp$Func))
combined_keggDF$Func_indx=as.numeric(combined_keggDF$Func)

ggplot()+
  geom_point(data = combined_keggDF, aes(as.factor(Func_indx), Trait, size=Enrichment, fill=-log10(Pvalue)), color="lightgrey", shape=21)+
  scale_fill_gradient(low="white", high="darkblue")+
  geom_point(data =combined_keggDF[combined_keggDF$`-log10(FDR)`>=1.3, ], aes(as.factor(Func_indx), Trait, size=Enrichment), shape=21, color="black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.5, colour = "black"))+
  theme(axis.text.y = element_blank(), panel.background = element_blank(), panel.grid = element_line(colour = "lightgrey"))+
  facet_wrap(~Parts, ncol = 1, scales = "free")+
  xlab("")+ylab("Traits")

tmpDF=unique(combined_keggDF[,c("Func_indx","Func")])
tmpDF$label=paste0(tmpDF$Func_indx," - ",tmpDF$Func)
tmpDF$hsa=gsub("\\s+.+","",tmpDF$Func) #%>% gsub("hsa","",.)
tmp_rbind=rbind(poolDF_new_ImmC_top100_funcs$dflist$KEGG,poolDF_new_EpiC_top100_funcs$dflist$KEGG,poolDF_new_StromC_top100_funcs$dflist$KEGG)
tmp_rbind=tmp_rbind[,c("ID","Description")] %>% unique()
colnames(tmp_rbind)=c("hsa","label")
tmpDF=left_join(tmpDF, tmp_rbind, by="hsa")
tmpDF$label=paste0(tmpDF$Func_indx," - ",tmpDF$label.y)
p=ggplot(tmpDF)+geom_text(aes(x=1, y=-Func_indx, label=label), hjust=0, color="black")+theme_nothing()
runFig(paste0(OutputDir,"/", paste0("Fig5.FuncEnrich_by3parts_legend_fullname")),p,8,20)
