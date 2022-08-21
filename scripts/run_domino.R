
library(Seurat)
library(domino)
library(optparse)

## reading parameters
option_list <- list(  
  make_option(c("-r", "--rds"), type="character", 
              help="rds (Seurat object) path"),
  make_option(c("-s", "--scenicdata"), type="character", 
              help="SCENIC output dir"),
  make_option(c("-o", "--output"), type="character",
              help="output dir"),
  make_option(c("-d", "--database"), type="character", default = '/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/test/cellphonedb_cellchat',
              help="database dir"),
  make_option(c("-t", "--cttrans"), type="character", default = NULL,
              help="celltype trans file"),
  make_option(c("-f", "--flag"), type="character", default = 'y',
              help="flag for Seurat object, whether it is in MAESTRO structure")
  
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

rds_path <- opts$rds
cpdb_dir <- opts$database
scenic_output_dir <- opts$scenicdata
domino_output_dir <- opts$output
ct_trans_file <- opts$cttrans
mastro_flag <- opts$flag


#### -------- start Domino workflow -------- ####

## loading data from Seurat object
RNA.res = readRDS(rds_path)
if (mastro_flag == 'y'){
  RNA.res = RNA.res$RNA
}
counts = as.matrix(RNA.res@assays$RNA@counts)
z_scores = RNA.res@assays$RNA@scale.data
clusters = as.factor(RNA.res@meta.data$assign.curated)
##### the next step only for human_heart or other datasets need to change celltype name ####
if (length(ct_trans_file) > 0){
  ct_curated = RNA.res@meta.data$assign.curated
  ct_trans_df = read.csv(ct_trans_file, sep='\t', check.names=F)
  ct_curated = plyr::mapvalues(x= ct_curated, from= as.vector(ct_trans_df[,2]), to=as.vector(ct_trans_df[,1]))
  clusters = as.factor(ct_curated)
}
#########################################################################

auc = t(read.table(paste0(scenic_output_dir,'/auc_mtx.csv'), header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE, sep = ','))


## create domino object
dom.obj = create_domino(signaling_db = cpdb_dir, 
                        features = auc, counts = counts, z_scores = z_scores, clusters = clusters, 
                        df = paste0(scenic_output_dir,'/regulons.csv'))


## run analysis using its default parameter (we set a larger max_tf_per_clust than its default, which default (in its example) is 10)
dom.res = build_domino(dom.obj, max_tf_per_clust = 20,
                       min_tf_pval = .001, max_rec_per_tf = 10, rec_tf_cor_threshold = .25)



########### -------------- extract results -------------- #########

## extract lr result matrix
signal_matrix_list = dom.res@cl_signaling_matrices
rec_ct_list = names(signal_matrix_list)

if (!file.exists(paste0(domino_output_dir,'/lr_ct_score'))){
  dir.create(paste0(domino_output_dir,'/lr_ct_score'))
}

for (ct in rec_ct_list){
  
  # extract TF
  mat = signal_matrix_list[[ct]]
  tfs = dom.res@linkages$clust_tf[[ct]]
  
  # TF to receptor
  all_recs = c()
  all_tfs = c()
  for(tf in tfs){
    recs = dom.res@linkages$tf_rec[[tf]]
    all_recs = c(all_recs, recs)
    if(length(recs)){
      all_tfs = c(all_tfs, tf)
    }
  }
  all_recs = unique(all_recs)
  all_tfs = unique(all_tfs)
  
  # Recs to ligs
  all_lrs = c()
  allowed_ligs = rownames(dom.res@cl_signaling_matrices[[ct]])
  for(rec in all_recs){
    tmp_ligs = dom.res@linkages$rec_lig[[rec]]
    common_lig = intersect(tmp_ligs, allowed_ligs)
    if (length(common_lig) > 0){
      for (c_lig in common_lig){
        lr_pair = paste0(c_lig, ' - ', rec)
        if (!lr_pair %in% all_lrs){
          if (length(all_lrs) == 0){
            # empty, init lr_ct_score_df (result matrix)
            lr_ct_score_df = mat[c_lig,]
          } else {
            lr_ct_score_df = rbind(lr_ct_score_df, mat[c_lig,])
          }
          all_lrs = c(all_lrs, lr_pair)
        }
      }
    }
    
  }
  
  rownames(lr_ct_score_df) = all_lrs
  colnames(lr_ct_score_df) = paste(gsub('L_','',colnames(mat)),ct,sep='|')
  
  write.table(lr_ct_score_df, file=paste0(domino_output_dir,'/lr_ct_score/',ct,'.tsv'), sep='\t', quote = F)
  
}


## save domino object
saveRDS(dom.res, file=paste(domino_output_dir,'dom_res.rds', sep = '/'))
saveRDS(dom.obj, file=paste(domino_output_dir,'dom_obj.rds', sep = '/'))



