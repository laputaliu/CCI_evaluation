

library(optparse)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(pheatmap)

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


run_cc <- function(count_path, meta_path, output_path){
  
  print('############ ------------- cellchat --------------- ############')
  
  print(paste0('>>> loading library and data <<< [', Sys.time(),']'))
  
  ######### prepare ##########
  data.norm = read.table(count_path,sep='\t',header=T,stringsAsFactors = F,row.names = 1)
  meta = read.table(meta_path,sep = '\t',header = F,stringsAsFactors = F,row.names = 1) 
  barcode = rownames(meta)
  label = as.vector(meta$V2)
  names(label) = barcode
  
  print(paste0('>>> start CellChat workflow <<< [', Sys.time(),']'))
  
  ########## start ######
  
  identity = data.frame(group = label, row.names = names(label)) # create a dataframe consisting of the cell labels
  unique(identity$group) # check the cell labels
  
  # Create a CellChat object
  cellchat <- createCellChat(object = as.matrix(data.norm), group.by='group', meta = identity)
  
  
  # Add cell information into meta slot of the object
  cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  # Preprocessing the expression data for cell-cell communication analysis
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multiprocess", workers = 1) # do parallel, but we can not do it on our own server
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  print(paste0('>>> Infer CCI network <<< [', Sys.time(),']'))
  
  
  ################ --------- Inference of cell-cell communication network -------------- ###################
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  
  # Infer the cell-cell communication at a signaling pathway leve
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  print(paste0('>>> saving results <<< [', Sys.time(),']'))
  
  if (substr(output_path, nchar(output_path), nchar(output_path)) == '/'){
    rda_output_path = paste0(output_path,'cellchat.rda')
  } else {
    rda_output_path = paste(output_path,'cellchat.rda',collapse = '/')
    output_path = paste0(output_path,'/')
  }
  
  save(cellchat,file = rda_output_path)
  
  
  ### output cellchat results (pval)
  cellchat_pval = cellchat@net$pval
  interaction_name = dimnames(cellchat_pval)[[3]]
  inter_types = c()
  pval_matrix = data.frame()
  for (i in 1:length(dimnames(cellchat_pval)[[1]])){
    for (j in 1:length(dimnames(cellchat_pval)[[2]])){
      inter_type = paste(dimnames(cellchat_pval)[[1]][i],dimnames(cellchat_pval)[[1]][j],sep = '|')
      inter_types = c(inter_types,inter_type)
      if (i+j == 2){
        pval_matrix = cellchat_pval[i,j,]
      }else{
        pval_matrix = cbind(pval_matrix,cellchat_pval[i,j,])
      }
    }
  }
  colnames(pval_matrix) = inter_types
  
  write.table(pval_matrix,file = paste0(output_path,'pval_matrix.tsv'), sep = '\t', quote = F)
  
  ### output cellchat results (prob)
  cellchat_prob = cellchat@net$prob
  interaction_name = dimnames(cellchat_prob)[[3]]
  inter_types = c()
  prob_matrix = data.frame()
  for (i in 1:length(dimnames(cellchat_prob)[[1]])){
    for (j in 1:length(dimnames(cellchat_prob)[[2]])){
      inter_type = paste(dimnames(cellchat_prob)[[1]][i],dimnames(cellchat_prob)[[1]][j],sep = '|')
      inter_types = c(inter_types,inter_type)
      if (i+j == 2){
        prob_matrix = cellchat_prob[i,j,]
      }else{
        prob_matrix = cbind(prob_matrix,cellchat_prob[i,j,])
      }
    }
  }
  colnames(prob_matrix) = inter_types
  
  write.table(prob_matrix,file = paste0(output_path, 'prob_matrix.tsv'), sep = '\t', quote = F)
  
  
  # output cellchat interaction names
  
  interaction_trans = cellchat@DB$interaction[,c('interaction_name','interaction_name_2')]
  write.table(interaction_trans, file = paste0(output_path, 'interaction_name_trans.tsv'),sep = '\t', quote = F, row.names = F)
  
  
  print(paste0('>>> plot heatmap <<< [', Sys.time(),']'))
  
  
  ################ ---------- count heatmap ------------ ################

  interaction_count = cellchat@net$count
  
  # rownames(cellchat@net$count)
  # colnames(interaction_count)
  # rownames(interaction_count)
  
  setwd(output_path)
  
  col1 = "dodgerblue4"
  col2 = 'peachpuff'
  col3 = 'deeppink4'
  col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
  show_rownames = T
  show_colnames = T
  scale="none"
  cluster_cols = T
  border_color='white'
  cluster_rows = T
  fontsize_row=11
  fontsize_col = 11
  main = ''
  treeheight_row=0
  family='Arial'
  treeheight_col = 0
  pheatmap(interaction_count, show_rownames = show_rownames, show_colnames = show_colnames, cluster_rows = F, cluster_cols = F,
           border_color=border_color, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
           main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = 'interaction_counts.pdf',display_numbers = F,fontsize_number = 11,number_color = 'white')
  
  pheatmap(log(interaction_count+1), show_rownames = show_rownames, show_colnames = show_colnames, cluster_rows = F, cluster_cols = F,
           border_color=border_color, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
           main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = 'interaction_counts_log.pdf',display_numbers = F,fontsize_number = 11,number_color = 'white')
  
  
  print(paste0('>>> done <<< [', Sys.time(),']'))
  
}


run_cc(count_path, meta_path, output_path)




