import pandas as pd
import os
import sys
import getopt


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -i input 
    -o output
    -c celltype
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

optlist,alist=getopt.getopt(sys.argv[1:],'hi:o:c:')
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-i':cytotalk_input_file = opt[1]
    elif opt[0] == '-o':output_dir = opt[1]
    elif opt[0] == '-c':done_ct_file = opt[1]

        
if __name__ == '__main__':

#     output_dir = '/fs/home/liuzhaoyang/project/cci_evaluation/human_heart/tools/CytoTalk/output'
#     cytotalk_input_file = '/fs/home/liuzhaoyang/project/cci_evaluation/human_heart/tools/CytoTalk/script/cytotalk.input'

    with open(done_ct_file,'r') as f:
        done_ct_list = [line.strip() for line in f.readlines()]
#     done_ct_list = [file for file in os.listdir(output_dir) if not file.startswith('.')]

    ct_input_df = pd.read_csv(cytotalk_input_file, sep='\t', header=None)
    ct_input_df.columns = ['cta','ctb','input','output','nohup']
    ct_input_df['pair'] = ct_input_df.apply(lambda x: '{}_{}'.format(x[0],x[1]),axis=1)

    select_index = [pair in done_ct_list for pair in list(ct_input_df['pair'])]
    reg_input_df = ct_input_df.loc[select_index,:]

    reg_input_out = '{}/cytotalk.input.acutal'.format('/'.join(cytotalk_input_file.split('/')[:-1]))
    reg_script_outdir = '{}/reg_part_script'.format('/'.join(cytotalk_input_file.split('/')[:-1]))

    if not os.path.exists(reg_script_outdir):
            os.mkdir(reg_script_outdir)

    script_name_list = []

    for i in range(reg_input_df.shape[0]):

        ct_a = reg_input_df.iloc[i,0]
        ct_b = reg_input_df.iloc[i,1]
        input_path = reg_input_df.iloc[i,2]
        output_path = reg_input_df.iloc[i,3]
        pair = reg_input_df.iloc[i,-1]

        part_script_name = '{}/regenerate_{}.R'.format(reg_script_outdir,pair)

        reg_script = '''rm(list = ls())

BetaUpperLimit <- 100
Order <- 5
no_cores <- 6

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


CellTypeA = '{}'
CellTypeB = '{}'
InputPath = '{}'
OutputPath = '{}'

print(paste0('>>>>>>>>>  ',CellTypeA, '|', CellTypeB,' <<<<<<<<<<'))

OutputPath_ct = paste0(OutputPath,CellTypeA,'_',CellTypeB)
if (file.exists(OutputPath_ct)){{

  #6) Generate the final signaling network between the two cell types. ~25min
  print(Sys.time())
  genSignalingNetwork(BetaUpperLimit, InputPath, OutputPath_ct,no_cores)

  #7) Extract ligand-receptor-associated pathways from the predicted signaling network between the two cell types.
  print(Sys.time())
  extractPathway(OutputPath_ct, Order)
  print(Sys.time())
  print("ALL DONE! CytoTalk has generated output files in the folders: 'OutputPath_ct/IllustratePCSF/' and 'OutputPath_ct/IllustratePathway/'")

  print(paste0('>>> end  <<< [', Sys.time(),']'))

    }}'''.format(ct_a, ct_b, input_path, output_path)
    
        with open(part_script_name,'w') as f:
            f.write(reg_script)

        script_name_list.append('regenerate_{}.R'.format(pair))


    with open('{}/cytotalk.part.script'.format('/'.join(cytotalk_input_file.split('/')[:-1])),'w') as f:
        f.write('\n'.join(script_name_list))


