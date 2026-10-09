# ============================================================================
# scibetRun helper and loading / merging of the public datasets
# ----------------------------------------------------------------------------
# Repository : 2026_SHCA_project
# Methods    : "Public scRNA-seq data processing"
# Source     : extracted from Analysis240525.r (L13173-13248, L13249-13267)
# Note       : analysis record. Paths point to the original data layout;
#              see the top-level README for the expected inputs.
# ============================================================================
scibetRun=function(trainfile,gene_info,outdir_prefix,celltype="major",nSelectGene=4000, pred_data=NULL){

  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment)
  library(data.table)
  library(tidyverse)
  library(tibble)
  library(scibet)

  gene_info = read.table("DATA_DIR",header = T)
  gene_info_sub = gene_info[gene_info$geneType %in% c('lncRNA',  'protein_coding','TR_C_gene','IG_C_gene', 'IG_C_pseudogene'),]
  f_train = trainfile
  train <- readRDS(f_train)
  train_set<- as_tibble(as(Matrix::t(train@assays$RNA$data), "dgCMatrix"))
  train_set['label'] = train@meta.data[,celltype]
  print(head(train_set))

  ############## Train model #################################
  genes = SelectGene(train_set, k = nSelectGene)
  genes = genes[genes %in% gene_info_sub$geneSymbol]
  gindex= sapply(genes,FUN=function(x){which(colnames(train_set)==x)})

  model = Learn(train_set, geneset = gindex-1)
  saveRDS(model, paste0(outdir_prefix,'.model.rds'))

  ############## Predict #################################
  read_data2list=function(infile, olist=NULL){

    if(!is.null(olist)){
      pred_data_list=olist
    }else{
      pred_data_list=list()
    }

    if(grepl(".rds",infile)){
      pred_data_list[[infile]]=readRDS(infile)
    }else if(grepl(".h5ad",infile)){
      sce <- readH5AD(infile, use_hdf5 = FALSE, reader = 'R')
      sce <- as.Seurat(sce, counts = 'X', data='X')
      sce@assays$RNA <- sce@assays$originalexp
      sce@assays$RNA$data <- log1p(sce@assays$RNA$data/100)
      pred_data_list[[infile]]=sce
    }

    return(pred_data_list)
  }

  pred_data_list=list()
  if(!is.null(pred_data)){
    if(typeof(pred_data)=="character"){
      if(grepl(".rds|.h5ad",pred_data)){
        pred_data_list=read_data2list(pred_data, pred_data_list)
      }else{
        pred_data="DATA_DIR"
        fff=list.files(pred_data,full.names = T,pattern = 'rds|h5ad')
        if(length(fff)==0){
          stop("unrecognized input; rds/h5ad or directory containing rds/h5ad supported!")
        }
        for(infile in fff[1:3]){
          pred_data_list=read_data2list(infile, pred_data_list)
        }
      }
    }

  }

}

########################## public dataset 202505 ###################
# NOTE: the original session sourced a personal package-loading helper here
#       (library(...) calls for Seurat, ggplot2, ComplexHeatmap, pheatmap, etc.).
#       Load the packages listed in the repository README instead.
setwd("./data/")

mypattern="Immune"
select_data_files=list.files(path = "DATA_DIR", pattern =  mypattern, full.names = T)

pub_data_list=list()
for(x in 1:length(select_data_files)){
pub_data_list[[x]]=readRDS(select_data_files[x])
}
pub_data=merge(pub_data_list[[1]], pub_data_list[2:length(pub_data_list)])
