import sys
import os
import getopt


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -c count matrix path
    -r raw count matrix path
    -m meta path
    -o tool output dir
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

optlist,alist=getopt.getopt(sys.argv[1:],'hc:r:m:o:')
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-c':count_path = opt[1]
    elif opt[0] == '-r':count_raw_path = opt[1]
    elif opt[0] == '-m':meta_path = opt[1]
    elif opt[0] == '-o':tool_output_dir = opt[1]


def generate_dir_script(count_path, count_raw_path, meta_path, tool_output_dir):
    '''
    generate the dir and scripts for each tool, (except iTALK)
    - dir : tool - {output, script}
    - script : named as submit_{tool_names}.sh, to run each tool, only need to run ./submit_{tool_names}.sh (maybe)
    
    the screen stuff can not accept the string that is too long, so if the code doesn't work, try to reduce the length of the string in stuff ""
    
    '''
    
    if not os.path.exists(tool_output_dir):
        os.mkdir(tool_output_dir)
    
    if tool_output_dir[-1] == '/':
        tool_output_dir = tool_output_dir[-1]

    ## prepare output dir

    tool_list = ['cc','cpdb','iTALK','SCR','NATMI','icellnet','NicheNet','CytoTalk','scmlnet','Connectome']

    for tool in tool_list:
        try:
            os.mkdir(os.path.join(tool_output_dir, tool))
            os.mkdir(os.path.join(tool_output_dir, tool, 'script'))
            os.mkdir(os.path.join(tool_output_dir, tool, 'output'))
        except FileExistsError:
            continue

      
    ## cellchat

    cc_command = 'screen -x cellchat -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_cc.R -c {} -m {} -o {}/cc/output/ > {}/cc/script/cc_nohup.out 2>&1 &\\n"'.format(count_path, meta_path, tool_output_dir, tool_output_dir)

    with open('{}/cc/script/submit_cc.sh'.format(tool_output_dir),'w') as f:
        f.write(cc_command)
        
    os.system('chmod 777 {}/cc/script/submit_cc.sh'.format(tool_output_dir))


    ## cpdb

    cpdb_command = 'screen -x cpdb -p 0 -X stuff "nohup time -v cellphonedb method statistical_analysis {} {} --counts-data gene_name --output-path {}/cpdb/output --threads 4 --database /fs/home/liuzhaoyang/project/cci_evaluation/Seurat_mapping/10X_mouse_brain/tools/cpdb/cellphonedb_user_2021-01-12-18_19.db > {}/cpdb/script/cpdb_nohup.out 2>&1 &\\n"'.format(tool_output_dir, meta_path, count_path, tool_output_dir, tool_output_dir)

    with open('{}/cpdb/script/submit_cpdb.sh'.format(tool_output_dir),'w') as f:
        f.write(cpdb_command)
        
    os.system('chmod 777 {}/cpdb/script/submit_cpdb.sh'.format(tool_output_dir))


    ## iTALK 直接在本地运行吧

    ## SCR

    scr_command = 'screen -x scr -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scr_cutoff_refine.R -c {} -m {} -o {}/SCR/output/ -s 0.6 > {}/SCR/script/scr_nohup.out 2>&1 &\\n"'.format(count_path, meta_path, tool_output_dir, tool_output_dir)

    with open('{}/SCR/script/submit_scr.sh'.format(tool_output_dir),'w') as f:
        f.write(scr_command)
        
    os.system('chmod 777 {}/SCR/script/submit_scr.sh'.format(tool_output_dir))



    ## NATMI

    natmi_command1 = 'sed \'1 iCell\\tAnnotation\' {} > {}/NATMI/meta_natmi.tsv && cd ~/biosoft/NATMI && time -v python ExtractEdges.py --interDB LRdb_cellchat --interSpecies human --emFile {} --annFile {}/NATMI/meta_natmi.tsv --species human --idType symbol --out {}/NATMI/output/'.format(meta_path, tool_output_dir, count_path, tool_output_dir, tool_output_dir)

    with open('{}/NATMI/script/run_natmi.sh'.format(tool_output_dir),'w') as f:
        f.write(natmi_command1)

    natmi_command = 'chmod 777 {}/NATMI/script/run_natmi.sh && screen -x natmi -p 0 -X stuff "cd {}/NATMI/script/ && nohup ./run_natmi.sh > {}/NATMI/script/natmi_nohup.out 2>&1 &\\n"'.format(tool_output_dir, tool_output_dir, tool_output_dir)

    with open('{}/NATMI/script/submit_natmi.sh'.format(tool_output_dir),'w') as f:
        f.write(natmi_command)
        
    os.system('chmod 777 {}/NATMI/script/submit_natmi.sh'.format(tool_output_dir))




    ## icellnet

    icellnet_command = 'screen -x icellnet -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_icellnet_refine.R -c {} -m {} -o {}/icellnet/output/ > {}/icellnet/script/icellnet_nohup.out 2>&1 &\\n"'.format(count_path, meta_path, tool_output_dir, tool_output_dir)

    with open('{}/icellnet/script/submit_icellnet.sh'.format(tool_output_dir),'w') as f:
        f.write(icellnet_command)
        
    os.system('chmod 777 {}/icellnet/script/submit_icellnet.sh'.format(tool_output_dir))


    ## NicheNet
    nichenet_command = 'screen -x nichenet -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_nichenet.R -c {} -m {} -o {}/NicheNet/output/ > {}/NicheNet/script/nichenet_nohup.out 2>&1 &\\n"'.format(count_path, meta_path, tool_output_dir, tool_output_dir)

    with open('{}/NicheNet/script/submit_nichenet.sh'.format(tool_output_dir),'w') as f:
        f.write(nichenet_command)
        
    os.system('chmod 777 {}/NicheNet/script/submit_nichenet.sh'.format(tool_output_dir))


    ## CytoTalk

    ## the prepare_input_cytotalk.py will generate all the files of the input dir and also for the run_cytotalk.sh and cytotalk.input
    cytotalk_command = 'python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_input_cytotalk.py -c {} -m {} -o {} && cd {}/CytoTalk/script/ && nohup time -v ./run_cytotalk.sh > {}/CytoTalk/script/cytotalk_nohup.out 2>&1 &'.format(count_path, meta_path, tool_output_dir, tool_output_dir, tool_output_dir)

    ## for the real run, using the submit_cytotalk.sh
    with open('{}/CytoTalk/script/submit_cytotalk.sh'.format(tool_output_dir),'w') as f:
        f.write(cytotalk_command)
        
    os.system('chmod 777 {}/CytoTalk/script/submit_cytotalk.sh'.format(tool_output_dir))
    
    ## if it occurs the error when generate the output, then we can use the script/regenerate_cytotalk_result.R



    ## scmlnet

    ## !!!!!!!!!!!!! scmlnet needs raw count data as input !!!!!!!!!!
    # scmlnet_command = 'screen -x scmlnet -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scmlnet.R -c {} -m {} -o {}/scmlnet/output/ -n {} > {}/scmlnet/script/scmlnet_nohup.out 2>&1 &\\n"'.format(count_raw_path, meta_path, tool_output_dir, 4, tool_output_dir)

    # with open('{}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir),'w') as f:
    #     f.write(scmlnet_command)

    # os.system('chmod 777 {}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir))

    scmlnet_command = 'sed \'1 iBarcode\tCluster\' {} > {}/scmlnet/scmlnet_meta.tsv && python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_input_scmlnet.py -c {} -m {}/scmlnet/scmlnet_meta.tsv -o {} -t 1 && cd {}/scmlnet/script/ && nohup time -v ./run_scmlnet.sh > {}/scmlnet/script/scmlnet_nohup.out 2>&1 &'.format(meta_path, tool_output_dir, count_raw_path, tool_output_dir, tool_output_dir, tool_output_dir, tool_output_dir)

    with open('{}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir),'w') as f:
        f.write(scmlnet_command)

    os.system('chmod 777 {}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir))




    ## Connectome
    connectome_command = 'screen -x connectome -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_connectome.R -c {} -m {} -o {}/Connectome/output/ > {}/Connectome/script/connectome_nohup.out 2>&1 &\\n"'.format(count_path, meta_path, tool_output_dir, tool_output_dir)

    with open('{}/Connectome/script/submit_connectome.sh'.format(tool_output_dir),'w') as f:
        f.write(connectome_command)
        
    os.system('chmod 777 {}/Connectome/script/submit_connectome.sh'.format(tool_output_dir))




if __name__ == '__main__':
    
#     print('{}\n{}\n{}\n{}\n'.format(count_path, count_raw_path, meta_path, tool_output_dir))
       
    generate_dir_script(count_path, count_raw_path, meta_path, tool_output_dir)


