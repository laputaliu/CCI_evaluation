

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
no_cores <- 4

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



run_cytotalk <- function(InputPath,OutputPath,WorkingPath,CellTypeA,CellTypeB,species,GeneFilterCutoff_A,GeneFilterCutoff_B,BetaUpperLimit,Order,no_cores){
  print(paste0('>>>>>>>>>  ',CellTypeA, '|', CellTypeB,' <<<<<<<<<<'))

  OutputPath_ct = paste0(OutputPath,CellTypeA,'_',CellTypeB)
  if (file.exists(OutputPath_ct)){
    
    result <- tryCatch({
      
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
  
}


run_cytotalk(InputPath,OutputPath,WorkingPath,CellTypeA,CellTypeB,species,GeneFilterCutoff_A,GeneFilterCutoff_B,BetaUpperLimit,Order,no_cores)



