
library(optparse)
library(Connectome)
library(ggplot2)
library(cowplot)
library(Seurat)


options(stringsAsFactors = FALSE)
option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count matrix / normalized count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltype annotation) path"),
  make_option(c("-o", "--output"), type="character",
              help="the output dir")
  )
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_path <- opts$count
meta_path <- opts$meta
output_path <- opts$output

if (!file.exists(output_path)){
  dir.create(output_path)
}


run_connectome <- function(count_path, meta_path, output_path){
  
  print('############ ------------- connectome --------------- ############')
  print(paste0('>>> loading data <<< [', Sys.time(),']'))
  count_df = read.csv(count_path, sep='\t', row.names=1)
  meta_df = read.table(meta_path,sep='\t')
  
  print(paste0('>>> generate seurat object <<< [', Sys.time(),']'))
  seuratObj <- CreateSeuratObject(counts = count_df)
  seuratObj@meta.data$celltype = meta_df$V2
  Idents(object = seuratObj) = meta_df$V2
  ct_list = levels(factor(meta_df$V2))

  print(paste0('>>> start Connectome workflow (preprocessing) <<< [', Sys.time(),']'))
  ## scale
  panc8 <- seuratObj
  connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
  genes <- connectome.genes[connectome.genes %in% rownames(panc8)]
  panc8 <- ScaleData(panc8,features = genes)
  
  print(paste0('>>> create Connectome object <<< [', Sys.time(),']'))
  panc8.con <- CreateConnectome(panc8,species = 'human',)
  
  ## interaction (edge) filtering
  print(paste0('>>> filter Connectome interactions <<< [', Sys.time(),']'))
  panc8.con2 <- FilterConnectome(panc8.con,min.pct = 0.1,min.z = 0.25,remove.na = T)
  
  print(paste0('>>> plot & save results <<< [', Sys.time(),']'))
  save(panc8.con2,file=file.path(output_path, 'CCI_network_filtered.rda'))
  output_df = panc8.con2$weight_sc
  names(output_df) = panc8.con2$edge
  write.table(output_df, file=file.path(output_path, 'weight_sc.tsv'),sep='\t',quote = F, col.names = F)
  print(paste0('>>> done <<< [', Sys.time(),']'))
  
}


run_connectome(count_path, meta_path, output_path)

