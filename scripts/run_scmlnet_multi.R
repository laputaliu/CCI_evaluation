
library(optparse)
library(Seurat)
library(Matrix)
library(parallel)
library(scMLnet)


options(stringsAsFactors = FALSE)

option_list <- list(  
  make_option(c("-a", "--ct_a"), type = "character",
              help = "LigClu"),
  make_option(c("-b", "--ct_b"), type = "character",
              help = "RecClu"),
  make_option(c("-c", "--count"), type="character", 
              help="count matrix / normalized count matrix path"),
  make_option(c("-m", "--meta"), type="character",
              help="meta data (celltypes annotation) path"),
  make_option(c("-o", "--output"), type="character",
              help="the output dir"),
  make_option(c("-n", "--ncores"), type = "integer",
              help = "the number of cores")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

LigClu <- opts$ct_a
RecClu <- opts$ct_b
count_path <- opts$count
meta_path <- opts$meta
output_path <- opts$output
ncores <- opts$ncores

if (!file.exists(output_path)){
  dir.create(output_path)
}

###### RunMLnet function #######
#### the default core number in scMLnet's original codes was set as the (max core number - 1). 
#### We rewrite some functions to make it possible to assign exact core number to avoid occupying all the computing resources while running on our public server.

getHighExpGene_new <- function(GCMat,barCluTable,CluNum.1,CluNum.2,pval,logfc,cores)
{

  #Function
  getNormData <- function(data.use){

    LogNorm <- function(data, scale_factor, display_progress = TRUE) {
      .Call('_Seurat_LogNorm', PACKAGE = 'Seurat', data, scale_factor, display_progress)
    }

    LogNormalize <- function(data, scale.factor = 1e4, display.progress = TRUE)
    {
      if (class(x = data) == "data.frame")
      {
        data <- as.matrix(x = data)
      }
      if (class(x = data) != "dgCMatrix") {
        data <- as(object = data, Class = "dgCMatrix")
      }
      # call Rcpp function to normalize
      if (display.progress) {
        cat("Performing log-normalization\n", file = stderr())
      }
      norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = display.progress)
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      return(norm.data)
    }

    data.test <- LogNormalize(data.use)

    return(data.test)

  }

  ChooseGene <- function(data.use,cells.1,cells.2,logfc){3

    #Function
    calpct <- function(data.use,cells,thresh.min=0){

      if(length(genes.use) > 20000)
      {
        num=ceiling(length(genes.use)/5000)

        data.temp <- c()
        for(i in seq(1,num,1))
        {
          start=1+5000*(i-1)
          if(i == num)
          {
            end=length(genes.use)
          }
          else
          {
            end=5000*i
          }

          tempGeneUse=genes.use[start:end]

          data.tempX <- round(
            x = apply(
              X = data.use[tempGeneUse, cells, drop = F],
              MARGIN = 1,
              FUN = function(x) {
                return(sum(x > thresh.min) / length(x = x))
              }
            ),
            digits = 3
          )#round
          data.temp <- c(data.temp,data.tempX)
        }#for
      }else{
        data.temp <- round(
          x = apply(
            X = data.use[, cells, drop = F],
            MARGIN = 1,
            FUN = function(x) {
              return(sum(x > thresh.min) / length(x = x))
            }
          ),
          digits = 3
        )#round
      }

      return(data.temp)
    }

    genes.use <- rownames(data.use)
    print(paste0("gene:",length(genes.use)))

    #choose Gene based on percent expressed
    min.pct <- 0.05
    min.diff.pct <- -Inf

    ##calculate pct
    data.temp1 <- calpct(data.use,cells.1)
    data.temp2 <- calpct(data.use,cells.2)
    data.alpha <- cbind(data.temp1, data.temp2)
    colnames(x = data.alpha) <- c("pct.1","pct.2")

    ##max pct
    alpha.max <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.max) <- rownames(x = data.alpha)
    genes.use <- names(x = which(x = alpha.max > min.pct))

    ##diff pct
    alpha.diff <- alpha.max - apply(X = data.alpha, MARGIN = 1, FUN = min)
    genes.use <- names(
      x = which(x = alpha.max > min.pct & alpha.diff > min.diff.pct)
    )

    #choose Gene based on average difference
    logfc.threshold <- logfc
    print("logfc.threshold:")
    print(logfc.threshold)
    pseudocount.use <- 1

    ##diff expr
    data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
    data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
    total.diff <- (data.1 - data.2)
    genes.diff <- names(x = which(x = total.diff > logfc.threshold))

    #filter
    genes.use <- intersect(x = genes.use, y = genes.diff)

    #output final gene Table(gene,logfc)
    geneTable <- data.frame(gene = genes.use,logfc = total.diff[genes.use])

    return(geneTable)
  }

  do.ttest <- function(data.test,genes.use,cells.1,cells.2,cores){

    #Function
    getpval <- function(gene){
      x <- as.vector(data.test[gene, cells.1])
      y <- as.vector(data.test[gene, cells.2])
      result <- t.test(x, y)
      return(result$p.value)
    }

    if(cores == 1){

      print("do T-test")
      p_val <- lapply(genes.use,getpval)
      p_val <- unlist(p_val)

    }else{

      print("T-test in parallel")
      cl <- makeCluster(cores)
      clusterEvalQ(cl, library(stats))
      clusterEvalQ(cl, library(Matrix))
      clusterExport(cl, varlist=c("data.test", "cells.1", "cells.2"), envir=environment())
      p_val <- parLapply(cl, genes.use, getpval)
      p_val <- unlist(p_val)
      stopCluster(cl)

    }

    gc()
    return(p_val)
  }

  #Main
  stri <- paste0("get High Exp Gene in ",CluNum.1)
  print(stri)
  barListResult <- getBarList(CluNum.1,CluNum.2,barCluTable)
  cells.1 <- barListResult[[1]]
  cells.2 <- barListResult[[2]]
  stri <- paste0(CluNum.1,":",length(cells.1))
  print(stri)

  data.use <- GCMat
  ChoGeneTable <- ChooseGene(data.use,cells.1,cells.2,logfc)
  genes.use <- rownames(ChoGeneTable)

  #normalize
  data.test <- getNormData(GCMat)
  data.test <- data.test[genes.use,c(cells.1,cells.2)]

  #do t test
  p_val <- do.ttest(data.test,genes.use,cells.1,cells.2,cores)

  #get final Gene
  GenePval <- data.frame(Gene = genes.use,LogFC = ChoGeneTable$logfc ,Pval = p_val,row.names = genes.use)
  GenePval1 <- GenePval[which(GenePval$Pval < pval),]

  # HighExpGene <- GenePval1$Gene
  # HighExpGene <- as.vector(HighExpGene)
  # HighExpGeneNum <- length(HighExpGene)
  #
  # gc()
  # stri <- paste("find high gene num:",HighExpGeneNum,sep="")
  # print(stri)
  # print("-----------------------")
  # return(HighExpGene)
  return(GenePval1)
}

source('/fs/home/liuzhaoyang/biosoft/scMLnet-master/R/Run_scMLnet.R')

run_scmlnet <- function(LigClu, RecClu, count_path, meta_path, output_path, ncores){
  
  print('############ ------------- scmlnet --------------- ############')
  print(paste0('>>> loading library and data <<< [', Sys.time(),']'))
  # scMLnet needs raw count data as input
  count_df = read.csv(count_path,sep='\t',row.names=1,header=T,check.names = F)
  count_df = as(as.matrix(count_df),"dgCMatrix")
  meta_df = read.table(meta_path,sep='\t',header = T)
  ct_list = levels(factor(meta_df$Cluster))
  
  if (substr(output_path, nchar(output_path), nchar(output_path)) != '/'){
    output_path = paste0(output_path,'/')
  }
  
  
  print(paste0('>>> start workflow for each cell type <<< [', Sys.time(),']'))
#   pval <- 0.05
#   logfc <- 0.15
  pval <- 0.1
  logfc <- 0.1
  LigRecLib <- "/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/scMLnet/database/LigRec.txt"
  TFTarLib <- "/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/scMLnet/database/TFTargetGene.txt"
  RecTFLib <- "/fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/scMLnet/database/RecTF.txt"
  LigRecFile <- LigRecLib
  TFTableFile <- TFTarLib
  RecTFTableFile <- RecTFLib
  cores=ncores
  
  print(paste0('>>> start ', LigClu, '_', RecClu,' <<< [', Sys.time(),']'))
  netList <- RunMLnet(count_df, meta_path, RecClu, LigClu, pval, logfc, LigRecLib, TFTarLib, RecTFLib,cores)
  write.table(netList$LigRec,file=paste0(output_path, LigClu, '_',RecClu, '_LR.txt'),sep='\t',quote=F,row.names = F,col.names = F)

  #remove cell not in GCMat
  BarCluTable <- read.table(meta_path,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  allCell <- colnames(count_df)
  tableCell <- rownames(meta_path)
  print("check table cell:")
  BarCluTable <- BarCluTable[which(BarCluTable$Barcode %in% allCell),]
  print(dim(BarCluTable))

  RecClus <- getHighExpGene_new(count_df,BarCluTable,RecClu,LigClu,pval,logfc,cores)
  LigClus <- getHighExpGene_new(count_df,BarCluTable,LigClu,RecClu,pval,logfc,cores)

  write.table(RecClus, file=paste0(output_path, LigClu, '_',RecClu, '_RecClus.tsv'),sep='\t',quote=F)
  write.table(LigClus, file=paste0(output_path, LigClu, '_',RecClu, '_LigClus.tsv'),sep='\t',quote=F)
  print(paste0('>>> end ', LigClu, '_', RecClu,' <<< [', Sys.time(),']')) 
  print(paste0('>>> all finished <<< [', Sys.time(),']'))
  
}


run_scmlnet(LigClu, RecClu, count_path, meta_path, output_path, ncores)

                    
