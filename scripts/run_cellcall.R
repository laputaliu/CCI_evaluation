
library(cellcall)
library(optparse)


## reading parameters
option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count file path"),
  make_option(c("-m", "--meta"), type="character", 
              help="meta (celltype) file path"),
  make_option(c("-o", "--output"), type="character",
              help="output dir"),
  make_option(c("-d", "--delim"), type="character", default = '@',
              help="the delimiter used to separate barcode and celltype in cell label"),
  make_option(c("-f", "--field"), type="integer", default = 2,
              help="the field of celltype in the cell label")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_file <- opts$count
meta_file <- opts$meta
output_dir <- opts$output
names.delim <- opts$delim
names.field <- opts$field


#### -------- start CellCall workflow -------- ####

### encode the celltype to the cell barcode
count_df = read.table(count_file, sep='\t', row.names = 1, header=T, check.names = F)

meta_df = read.table(meta_file, sep='\t', row.names = 1, check.names = F)
meta_df = meta_df[colnames(count_df),]
# since the cell id cannot have '-/_', so we replace '-' using '.', '_' using '..', and we will change them back in the final result
meta_df = gsub('-','.',meta_df)
meta_df = gsub('_','..',meta_df)
update_barcode = paste(colnames(count_df),meta_df, sep=names.delim)
colnames(count_df) = update_barcode


### create cellcall object
cellcall.obj <- CreateNichConObject(data=count_df, min.feature = 0,
                                    names.field = names.field,
                                    names.delim = names.delim,
                                    source = "UMI",
                                    scale.factor = 10^6,
                                    Org = "Homo sapiens",
                                    project = "cellcall")


### infer cci
cellcall.obj <- TransCommuProfile(object = cellcall.obj,
                                  pValueCor = 0.05,
                                  CorValue = 0.1,
                                  topTargetCor=1,
                                  p.adjust = 0.05,
                                  use.type="median",
                                  probs = 0.9,
                                  method="weighted",
                                  IS_core = TRUE,
                                  Org = 'Homo sapiens')


### output result
# expr_l_r_log2_scale: The score of ligand-receptor in cellA-cellB with log transform and scale to [0,1]
lr_exp_log2_scale <- cellcall.obj@data$expr_l_r_log2_scale
lr_cols = colnames(lr_exp_log2_scale)
lr_cols = gsub('-','|',lr_cols)
lr_cols = gsub('\\.\\.','_',lr_cols)
lr_cols = gsub('\\.','-',lr_cols)
lr_rows = rownames(lr_exp_log2_scale)
lr_rows = gsub('-',' - ',lr_rows)
colnames(lr_exp_log2_scale) = lr_cols
rownames(lr_exp_log2_scale) = lr_rows

write.table(lr_exp_log2_scale, file=paste0(output_dir,'/','lr_exp_log2_scale.tsv'), sep='\t', quote = F)
saveRDS(cellcall.obj, file=paste0(output_dir,'/','cellcall.obj.rds'))





