


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
              help="meta data (celltypes annotation) path"),
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
  print(paste0('>>> load library and data <<< [', Sys.time(),']'))
  
  ## local path
  count_df = read.csv(count_path, sep='\t', row.names=1)
  meta_df = read.table(meta_path,sep='\t')
  
  if (substr(output_path, nchar(output_path), nchar(output_path)) != '/'){
    output_path = paste0(output_path,'/')
  }
  
  print(paste0('>>> generate seurat object <<< [', Sys.time(),']'))
  
  seuratObj <- CreateSeuratObject(counts = count_df)
  seuratObj@meta.data$celltype = meta_df$V2
  Idents(object = seuratObj) = meta_df$V2
  
  ct_list = levels(factor(meta_df$V2))

  print(paste0('>>> start Connectome workflow (preprocessing) <<< [', Sys.time(),']'))
  
  ##### scale
  panc8 <- seuratObj
  connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
  genes <- connectome.genes[connectome.genes %in% rownames(panc8)]
  panc8 <- ScaleData(panc8,features = genes)
  ## min.cells.per.ident set to 5 , the cell(spot) number of Endothelial is just 8
  
  print(paste0('>>> create Connectome object <<< [', Sys.time(),']'))
  
  panc8.con <- CreateConnectome(panc8,species = 'human',)
  
  print(paste0('>>> filter Connectome interactions <<< [', Sys.time(),']'))
  
  ## interaction (edge) filtering
  panc8.con2 <- FilterConnectome(panc8.con,min.pct = 0.1,min.z = 0.25,remove.na = T)
  # min.pct = 0.1, min.z = 0.5  -->  80 edges
  # min.pct = 0.1,min.z = 0.25  -->  720 edges
  
  print(paste0('>>> analysis finished <<< [', Sys.time(),']'))
  
  print(paste0('>>> plot & save results <<< [', Sys.time(),']'))
  
  save(panc8.con2,file=paste0(output_path, 'CCI_network_filtered.rda'))
  
  output_df = panc8.con2$weight_sc
  names(output_df) = panc8.con2$edge
  write.table(output_df, file=paste0(output_path, 'weight_sc.tsv'),sep='\t',quote = F, col.names = F)
  
  print(paste0('>>> done <<< [', Sys.time(),']'))
  
}


run_connectome(count_path, meta_path, output_path)

















