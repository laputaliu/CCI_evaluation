
library(iTALK)
library(magrittr)
library(dplyr)
library(optparse)

options(stringsAsFactors = FALSE)
option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count matrix / normalized count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltypes annotation) path"),
  make_option(c("-d", "--database"), type="character",
              help="iTALK database path"),
  make_option(c("-n", "--ntop"), type="character",
              help="top_genes number", default = 50),
  make_option(c("-o", "--output"), type="character",
              help="the output dir")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_path <- opts$count
meta_path <- opts$meta
db_path <- opts$database
ntop <- opts$ntop
output_path <- opts$output

if (!file.exists(output_path)){
  dir.create(output_path)
}

run_italk <- function(count_path, meta_path, db_path, ntop, output_path){
  
  print('############ ------------- italk --------------- ############')
  print(paste0('>>> loading data <<< [', Sys.time(),']'))
  load(db_path)
  comm_list<-c('growth factor','other','cytokine','checkpoint')
  res<-NULL
  
  # trandform input count & meta data to required input format: row as cell, col as genes, the last column is cell_type
  data<-read.table(count_path, sep='\t', header=T, stringsAsFactors = F, row.names=1, check.names = F)
  data <- t(data)
  cell_type = read.table(meta_path, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, check.names = F) 
  colnames(cell_type) = 'cell_type'
  data = cbind(data,cell_type)
  
  setwd(output_path)
  
  print(paste0('>>> start italk workflow <<< [', Sys.time(),']'))
  highly_exprs_genes<-rawParse(data,top_genes=ntop,stats='mean')
  res<-FindLR(highly_exprs_genes,datatype='mean count',comm_type="other", database = database_cellchat)

  print(paste0('>>> write result <<< [', Sys.time(),']'))
  res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),]
  saveRDS(res, file = 'FindLR_P2_rep2_st_cellchatdb_res.rds')
  write.table(res,file='LR_result.tsv',sep='\t',quote = F,row.names = F)
  saveRDS(highly_exprs_genes,'highly_exprs_genes.rds')
  print(paste0('>>> all done <<< [', Sys.time(),']'))
}

run_italk(count_path, meta_path, db_path, ntop, output_path)

