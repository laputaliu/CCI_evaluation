library(Seurat)
library(optparse)

options(stringsAsFactors = FALSE)

option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltypes annotation) path"),
  make_option(c("-o", "--output"), type="character",
              help="the output path")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_path <- opts$count
meta_path <- opts$meta
output_path <- opts$output


count.data <- read.csv(file=count_path,sep='\t',row.names = 1,check.names = F)

obj <- CreateSeuratObject(counts = count.data, min.cells = 3, project = "Domino")
obj <- NormalizeData(obj)
obj <- ScaleData(obj, features=rownames(obj))

meta_df = read.csv(file=meta_path,sep='\t',row.names = 1,check.names = F, header=F)
ct_curated = as.character(meta_df$V2)

obj@meta.data$assign.curated = ct_curated

output_dir = unlist(strsplit(output_path,split = '/'))
output_dir = output_dir[-length(output_dir)]
output_dir = paste(output_dir, collapse ='/')

if (!file.exists(output_dir)){
    dir.create(output_dir,recursive=TRUE)
}

saveRDS(obj, file = output_path)