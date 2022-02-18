library(optparse)
library(SingleCellSignalR)

options(stringsAsFactors = FALSE)

option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count matrix / normalized count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltypes annotation) path"),
  make_option(c("-o", "--output"), type="character",
              help="the output dir"),
  make_option(c("-s", "--scorecutoff"), type="double",
              help="the score cutoff")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_path <- opts$count
meta_path <- opts$meta
output_path <- opts$output
s.score <- opts$scorecutoff

if (!file.exists(output_path)){
  dir.create(output_path)
}


run_scr <- function(count_path, meta_path, output_path,s.score=0.6){
  
  print('############ ------------- scr --------------- ############')
  print(paste0('>>> load library and data <<< [', Sys.time(),']'))
  
  if (substr(output_path, nchar(output_path), nchar(output_path)) != '/'){
    output_path = paste0(output_path,'/')
  }
  
  data = read.csv(count_path, sep = '\t',row.names=1)
  genes = rownames(data)
  
  ## celltype names to cluster number
  meta_df = read.table(meta_path, sep='\t')
  celltype_list = levels(as.factor(meta_df$V2))
  cluster = meta_df$V2
  for (n in seq(1,length(celltype_list))){
    ct = celltype_list[n]
    cluster = replace(cluster, cluster == ct, n)
  }
  cluster = as.numeric(cluster)
  
  # ###--- codes below are for generating scr's LRdb from cellchatdb, which only need to be run once. ---###
  # print(paste0('>>> generate LR db <<< [', Sys.time(),']'))
  # ##### generate cellchatdb #####
  # ligand = unlist(read.table('/fs/home/liuzhaoyang/project/cci_evaluation/SCR/database_cellchat/lignad_cellchat.tsv'))
  # receptor = unlist(read.table('/fs/home/liuzhaoyang/project/cci_evaluation/SCR/database_cellchat/receptor_cellchat.tsv'))
  # names(ligand) = NULL
  # names(receptor) = NULL
  # LRdb = data.frame("ligand" = ligand,
  #                   "receptor" = receptor,
  #                   "ligand.name" = ligand,
  #                   "receptor.name" = receptor,
  #                   "ligand.synonyms" = ligand,
  #                   "receptor.synonyms" = receptor,
  #                   "ligand.altern.names" = ligand,
  #                   "receptor.altern.names" = receptor,
  #                   "source" = rep("cellchat",length(ligand)),
  #                   "PMIDs" = rep("PMID",length(ligand)),
  #                   "cells.L" = rep("",length(ligand)),
  #                   "cells.R" = rep("",length(ligand)),
  #                   "remarks" = rep("",length(ligand))
  # )
  # save(LRdb,file = paste0(output_path, "LRdb.rda"))
  
  load('/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/SCR/database_cellchat/LRdb.rda')
  
  print(paste0('>>> start SCR workflow <<< [', Sys.time(),']'))
  celltype_list = gsub(' ', '_', celltype_list)
  celltype_list = gsub('/', '_', celltype_list)
  
  # different expression genes
  clust.ana <- cluster_analysis(data = data, genes = genes, cluster = cluster, c.names = celltype_list, write=FALSE)
  
  signal <- cell_signaling(data = data, genes = genes, cluster = cluster, species = "homo sapiens", c.names = celltype_list,int.type = c("autocrine"), s.score = s.score, write=FALSE)
  # default s.score cutoff -> 0.5
  
  print(paste0('>>> write result <<< [', Sys.time(),']'))
  for (n in seq(1,length(signal))){
    ct_type = paste0(colnames(signal[[n]])[1],'_', colnames(signal[[n]])[2])
    ip_df = signal[[n]]
    ct_type = gsub('/','_', ct_type)
    ct_type = gsub(' ','_',ct_type)
    write.table(ip_df, paste0(output_path,ct_type,'.tsv'),sep='\t', quote=F, row.names = F)
  }
  
  print(paste0('>>> all done <<< [', Sys.time(),']'))
  
}

run_scr(count_path, meta_path, output_path, s.score)
