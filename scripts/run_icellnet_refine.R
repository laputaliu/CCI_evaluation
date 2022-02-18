library(optparse)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
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


LR.selection <- function (lr = lr, thresh = 0, topn = NULL, sort.by = "sum", db.name.couple=db.name.couple) {
  interactions = rownames(lr)
  lr = as.data.frame(lr)
  lr$Pair = interactions
  lr = as.data.frame(lr[stats::complete.cases(lr), ])
  lr = lr %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. >
                                                             0)) %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. >=
                                                                                                                    thresh))
  if (!is.null(topn)) {
    if (topn > dim(lr)[1]) {
      note(paste0("lr contains only ", dim(lr)[1], " after filtering interaction highest than theshold"))
    }
  }
  if (sort.by == "sum") {
    lr = lr %>% dplyr::mutate(sum = rowSums(dplyr::across(where(is.numeric)))) %>%
      dplyr::arrange(dplyr::desc(sum)) %>% dplyr::top_n(topn,
                                                        sum)
    lr = lr %>% dplyr::select(-sum)
  }
  else if (sort.by == "var") {
    lr$variance = apply(dplyr::select_if(lr, is.numeric),
                        1, var, na.rm = TRUE)
    lr = lr %>% dplyr::arrange(dplyr::desc(variance)) %>%
      dplyr::top_n(topn, variance)
    lr = lr %>% dplyr::select(-variance)
  } else stop("sort.by argument should be fixed on var or sum")
  
  rownames(lr)=lr$Pair
  lr =lr %>% dplyr::select(-Pair)
  
  return(lr)
  
}

run_icellnet <- function(count_path, meta_path, output_path){
  
  print('############ ------------- icellnet --------------- ############')
  
  ### loading database (cellchatdb)
  print(paste0('>>> loading data <<< [', Sys.time(),']'))
  db=as.data.frame(read.csv('/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/ICELLNET/database/cellchatdb.tsv', sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
  count_df = read.csv(count_path,sep='\t',row.names = 1,check.names = F)
  meta_df = read.table(meta_path,sep='\t',check.names = F)
  
  ### generate input seurat object
  print(paste0('>>> generate Seurat object <<< [', Sys.time(),']'))
  seuratObj <- CreateSeuratObject(counts = count_df)
  seuratObj@meta.data$celltype = meta_df$V2
  Idents(object = seuratObj) = meta_df$V2
  
  
  print(paste0('>>> start ICELLNET workflow (sc.data.cleaning) <<< [', Sys.time(),']'))
  # filter.perc=10
  filter.perc = 0
  average.clean= sc.data.cleaning(object = seuratObj, db = db, filter.perc = filter.perc, save_file = T, path=output_path, force.file = F)
  ct_list = levels(factor(meta_df$V2))
  data.icell=as.data.frame(gene.scaling(as.data.frame(average.clean), n=1, db=db))
  
  ## go through each celltype pair
  print(paste0('>>> Go through each cell types <<< [', Sys.time(),']'))
  output_ct_list = c()
  for (i in seq(length(ct_list))){
    ct_a = ct_list[i]
    
    PC = ct_list[which(ct_list != ct_a)]
    PC.data=as.data.frame(data.icell[,c(PC,"Symbol")], row.names = rownames(data.icell))
    
    PC.target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), "Class"=colnames(data.icell))
    rownames(PC.target)=PC.target$ID
    
    score.computation.1= icellnet.score(direction="out", PC.data=PC.data, 
                                        CC.data= as.data.frame(data.icell[,c(ct_a)], row.names = rownames(data.icell)),  
                                        PC.target = PC.target, PC=PC, CC.type = "RNAseq", 
                                        PC.type = "RNAseq",  db = db)
    score1=as.data.frame(score.computation.1[[1]])
    lr1=score.computation.1[[2]]
    lr1 = lr1[!duplicated(rownames(lr1)),]
    
    ## select the top LR pairs
    pairs=LR.selection(lr = lr1, thresh = 0 , topn=200 , sort.by="var")
    print(paste0('>>> ', ct_a, 'vs others finished <<< [', Sys.time(),']'))
    # colnames(pairs) = PC
    
    out_ct_a = gsub(' ','_',ct_a)
    out_ct_a = gsub('/','_',out_ct_a)
    
    write.table(as.matrix(pairs), file=paste0(output_path,'score_',out_ct_a,'_out_all.tsv'), sep='\t', quote=F)
  }
}


run_icellnet(count_path, meta_path, output_path)


