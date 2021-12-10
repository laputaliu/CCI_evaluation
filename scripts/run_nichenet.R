
library(optparse)
library(nichenetr)
library(Seurat)
library(tidyverse)

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


run_nichenet <- function(count_path, meta_path, output_path){

  print('############ ------------- nichenet --------------- ############')
  
  print(paste0('>>> loading library and data <<< [', Sys.time(),']'))
  
  ###### server ######
  count_df = read.table(count_path,sep='\t')
  meta_df = read.table(meta_path,sep='\t', row.names = 1)

  ligand_target_matrix = readRDS('/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/NicheNet/PAAD/input/ligand_target_matrix.rds')
  lr_network = readRDS('/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/NicheNet/PAAD/input/lr_network.rds')
  weighted_networks = readRDS('/fs/home/liuzhaoyang/project/CCI_tools/cci_evaluation/NicheNet/PAAD/input/weighted_networks.rds')
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
  print(paste0('>>> generate seurat object <<< [', Sys.time(),']'))
  
  ## try to generate Seurat object
  seuratObj <- CreateSeuratObject(counts = count_df)
  seuratObj@meta.data$celltype = meta_df$V2
  Idents(object = seuratObj) = meta_df$V2
  
  ct_list = levels(factor(meta_df$V2))

  print(paste0('>>> start Nichenet workflow for each cell types <<< [', Sys.time(),']'))
  
  
  ###### using for roop to go over each ct one by one ########
  for (ct_a in ct_list){
    for (ct_b in ct_list){
      if (ct_a == ct_b){
        next
      }
      
      print(paste0('>>> ', ct_a, '_', ct_b, ' start <<< [', Sys.time(), ']'))
      
      out_ct_a = gsub(' ','_',ct_a)
      out_ct_a = gsub('/','_',out_ct_a)
      out_ct_b = gsub(' ','_',ct_b)
      out_ct_b = gsub('/','_',out_ct_b)
      
      OutputPath = paste0(output_path,'output',out_ct_a,'_',out_ct_b)
      dir.create(OutputPath)
      
      receiver = ct_a
      expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
      background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
      
      sender_celltypes = c(ct_b)
      list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
      
      seurat_obj_receiver= subset(seuratObj, idents = receiver)
      
      ### Define a set of potential ligands
      
      ligands = lr_network %>% pull(from) %>% unique()
      receptors = lr_network %>% pull(to) %>% unique()
      
      expressed_ligands = intersect(ligands,expressed_genes_sender)
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      
      potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
      
      
      ### NicheNet ligand activity analysis
      ## rank the potential ligands based on the presence of their target genes in the gene set of interest
      ligand_activities = predict_ligand_activities(expressed_genes_sender, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
      ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
      
      #### get top 20 activity ligands
      best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
      
      ## plot senders and top ligands expression situation
      # pdf(file = '/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/sender&top_ligand_dotplot.pdf')
      # pdf(file = paste0(OutputPath,'/sender&top_ligand_dotplot.pdf'))
      # DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
      # dev.off()
      
      
      ##### Infer receptors and top-predicted target genes of top-ranked ligands
      
      ### target gene
      active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = expressed_genes_sender, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
      
      active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
      
      order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
      order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
      rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
      colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
      
      vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
      
      # pdf(file='/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/toprank_ligand&target_heatmap.pdf')
      # pdf(file=paste0(OutputPath, '/toprank_ligand&target_heatmap.pdf'))
      # p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
      # p_ligand_target_network
      # dev.off()
      
      ### receptor gene
      lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
      best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
      
      lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
      
      lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
      lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
      
      dist_receptors = dist(lr_network_top_matrix, method = "binary")
      hclust_receptors = hclust(dist_receptors, method = "ward.D2")
      order_receptors = hclust_receptors$labels[hclust_receptors$order]
      
      dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
      hclust_ligands = hclust(dist_ligands, method = "ward.D2")
      order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
      
      order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
      order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
      
      vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
      rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
      colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
      
      # pdf(file='/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/toprank_ligand&receptor_heatmap.pdf')
      # pdf(file=paste0(OutputPath, '/toprank_ligand&receptor_heatmap.pdf'))
      # p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
      # p_ligand_receptor_network
      # dev.off()
      
      print(paste0('>>> ', ct_a, '_', ct_b, ' write <<< [', Sys.time(), ']'))
      # write.table(vis_ligand_receptor_network, file='/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/LR_potential.tsv',quote=F, sep='\t')
      write.table(vis_ligand_receptor_network, file=paste0(OutputPath, '/LR_potential.tsv'),quote=F, sep='\t')
      
      print(paste0('>>> ', ct_a, '_', ct_b, ' finish <<< [', Sys.time(), ']'))
      
      
      ###### more strict ligand-receptor result
      ##### considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
      
      lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
      ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
      receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

      lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
      lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

      lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
      lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

      dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
      hclust_receptors = hclust(dist_receptors, method = "ward.D2")
      order_receptors = hclust_receptors$labels[hclust_receptors$order]

      dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
      hclust_ligands = hclust(dist_ligands, method = "ward.D2")
      order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

      order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
      order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

      vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
      rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
      colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

      # pdf(file='/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/toprank_ligand&receptor_strict_heatmap.pdf')
      # pdf(file=paste0(OutputPath, '/toprank_ligand&receptor_strict_heatmap.pdf'))
      # p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
      # p_ligand_receptor_network_strict
      # dev.off()
      print(paste0('>>> ', ct_a, '_', ct_b, ' strict write <<< [', Sys.time(), ']'))
      # write.table(vis_ligand_receptor_network, file='/fs/home/liuzhaoyang/project/cci_evaluation/NicheNet/PAAD/LR_potential_strict.tsv',quote=F, sep='\t')
      write.table(vis_ligand_receptor_network, file=paste0(OutputPath, '/LR_potential_strict.tsv'),quote=F, sep='\t')
      print(paste0('>>> ', ct_a, '_', ct_b, ' strict finish <<< [', Sys.time(), ']'))
      
    } 
}

}


run_nichenet(count_path, meta_path, output_path)













