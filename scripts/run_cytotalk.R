

################ ----------------------------------------- #################

## you need to split your count matrix by cell types first
## we finish this step in python, the code of this step can be found in our markdown note of CytoTalk.

## also, don't forget to copy deter_W4.R, gen_PCSF.py, LigandReceptor_Human.txt, ProteinCodingGenes_Human.txt these 
##   four files into your input dir.

## then, for running Cytotalk, we only need to source this script in R. e.g. source('run_Cytotalk.R')

############################################################################





### The objective of this R script: infer a signaling network between two specified cell types using a step-by-step running fasion.
# System requirements: >= 4 logical cores, >= 8GB memory, >= 10GB disk.
# Timing: ~6h (2.3 GHz 8-Core Intel Core i9, 16 logical cores and 16GB memory)


# There are 10 input parameters for running the CytoTalk algorithm as below.
#' @param species indicating the species from which the scRNA-Seq data are generated. Currently, “Human” and “Mouse” are supported.
#' @param CellTypeA specifying the name of cell type A. Please make sure that the cell type name should be consistent with scRNA-seq data file.
#' @param CellTypeB specifying the name of cell type B. Please make sure that the cell type name should be consistent with scRNA-seq data file.
#' @param GeneFilterCutoff_A indicating the cutoff for removing lowly-expressed genes in the processing of scRNA-seq data of cell type A. The default cutoff value is 0.1, which means that genes expressed in less than 10% of all cells of cell type A are removed.
#' @param GeneFilterCutoff_B indicating the cutoff for removing lowly-expressed genes in the processing of scRNA-seq data of cell type B. The default cutoff value is 0.1, which means that genes expressed in less than 10% of all cells of cell type B are removed.
#' @param BetaUpperLimit indicating the upper limit of the test values of the PCSF objective function parameter β, which is inversely proportional to the total number of genes in a given cell-type pair after removing lowly-expressed genes. Based on preliminary tests, the upper limit of β value is suggested to be 100 (default) if the total number of genes in a given cell-type pair is above 10,000. However, if the total number of genes is below 5000, it is necessary to increase the upper limit of β value to 500.
#' @param Order indicating the order of the neighborhood of the ligand and receptor to extract their-associated pathway from the predicted signaling network between the two cell types.
#' @param InputPath specifying the directory that includes necessary input files for running CytoTalk.
#' @param OutputPath specifying the directory that includes the folder "/IllustratePCSF/" and "/IllustratePathway/" to show the whole predicted network and each ligand-receptor-associated pathway, respectively.
#' @param WorkingPath specifying the directory that includes all CytoTalk functions.



####-------- running CytoTalk -----------####

### Customize your inputs and parameters:
rm(list = ls())


library(optparse)

option_list <- list(  
  make_option(c("-a", "--celltypeA"), type="character", 
              help="CellTypeA"),
  make_option(c("-b", "--cellTypeB"), type="character",
              help="CellTypeB"),
  make_option(c("-i", "--inputDir"), type="character",
              help="input dir"),
  make_option(c("-o", "--outputDir"), type="character",
              help="output dir")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

CellTypeA <- opts$celltypeA
CellTypeB <- opts$cellTypeB
InputPath <- opts$inputDir
OutputPath <- opts$outputDir

species <- "Human"
GeneFilterCutoff_A <- 0.1
GeneFilterCutoff_B <- 0.1
BetaUpperLimit <- 100
Order <- 5
no_cores <- 6

# InputPath <- "/fs/home/liuzhaoyang/project/cci_evaluation/human_heart/tools/CytoTalk/input"
# OutputPath <- "/fs/home/liuzhaoyang/project/cci_evaluation/human_heart/tools/CytoTalk/output/"

WorkingPath <- "/fs/home/liuzhaoyang/biosoft/CytoTalk_package_v3.1.0/CytoTalk_Function"
setwd(WorkingPath)

#----------------Import functions (DO NOT change code below)-----------------#
source("gen_intracellularNetMat.R")
source("comp_MIcoexp_TypA_WinPara.R")
source("comp_MIcoexp_TypB_WinPara.R")
source("comp_GeneNet_TypA_LinuxPara.R")
source("comp_GeneNet_TypB_LinuxPara.R")
source("comp_NonSelfTalkScore_TypA.R")
source("comp_NonSelfTalkScore_TypB.R")
source("construct_integratedNetwork.R")
source("JobRun_Parallel.R")
source("gen_signalingNetwork.R")
source("gen_signalingPathway.R")
#----------------Import functions (DO NOT change code above)-----------------#


run_cytotalk <- function(InputPath,OutputPath,WorkingPath,CellTypeA,CellTypeB,species,GeneFilterCutoff_A,GeneFilterCutoff_B,BetaUpperLimit,Order,no_cores){
  print(paste0('>>>>>>>>>  ',CellTypeA, '|', CellTypeB,' <<<<<<<<<<'))
  print(paste0('>>> load library and data <<< [', Sys.time(),']'))

  print(paste0('>>> start CytoTalk workflow <<< [', Sys.time(),']'))
  
  if (!file.exists(OutputPath)){
    dir.create(OutputPath)
  }
  
  if (substr(OutputPath, nchar(OutputPath), nchar(OutputPath)) != '/'){
    OutputPath = paste0(OutputPath,'/')
  }
  
  OutputPath_ct = paste0(OutputPath,CellTypeA,'_',CellTypeB)
  if (!file.exists(OutputPath_ct)){
    dir.create(OutputPath_ct)
  }

  #1) Pre-process the scRNA-seq data. ~5min
  # print("Start running CytoTalk!")
  # print(Sys.time())
  preprocess(species, CellTypeA, CellTypeB, GeneFilterCutoff_A, GeneFilterCutoff_B, InputPath, OutputPath_ct)
  
  result <- tryCatch({
    #2) Generate mutual information matrix using parallel environment. ~4h depending on the number of single cells and available logical cores.
    print(Sys.time())
    compMI_TypeA(OutputPath_ct,no_cores)
    print(Sys.time())
    compMI_TypeB(OutputPath_ct,no_cores)
    
    #3) Generate indirect edge-filtered network matrix using parallel environment. ~10min
    print(Sys.time())
    compGeneNet_TypeA(OutputPath_ct,no_cores)
    print(Sys.time())
    compGeneNet_TypeB(OutputPath_ct,no_cores)
    
    #4) Construct the integrated gene network. ~40min
    print(Sys.time())
    compNonSelfTalkScore_TypeA(species, CellTypeA, InputPath, OutputPath_ct)
    compNonSelfTalkScore_TypeB(species, CellTypeB, InputPath, OutputPath_ct)
    constructIntegratedNetwork(species, BetaUpperLimit, InputPath, OutputPath_ct)
    
    #5) Generate background PCSFs based on the integrated gene network. ~20min
    print(Sys.time())
    JobRunParallel(OutputPath_ct,no_cores)
    
    #6) Generate the final signaling network between the two cell types. ~25min
    print(Sys.time())
    genSignalingNetwork(BetaUpperLimit, InputPath, OutputPath_ct,no_cores)
    
    #7) Extract ligand-receptor-associated pathways from the predicted signaling network between the two cell types.
    print(Sys.time())
    extractPathway(OutputPath_ct, Order)
    print(Sys.time())
    print("ALL DONE! CytoTalk has generated output files in the folders: 'OutputPath_ct/IllustratePCSF/' and 'OutputPath_ct/IllustratePathway/'")
    
    print(paste0('>>> end  <<< [', Sys.time(),']'))
  }, warning = function(w) {
    print(e)
  }, error   = function(e) { 
    print(e)
    print(paste0('EEEEEEEEEEEEEEEE ---- ',CellTypeA,'|',CellTypeB,' ---- EEEEEEEEEEEEEEEE'))
  }, finally = {
    print(paste0('>>> end  <<< [', Sys.time(),']'))
  })
  
}


run_cytotalk(InputPath,OutputPath,WorkingPath,CellTypeA,CellTypeB,species,GeneFilterCutoff_A,GeneFilterCutoff_B,BetaUpperLimit,Order,no_cores)




