import sys
import os
import getopt
import pandas as pd
import itertools


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -c count matrix path
    -m meta path
    -o tool output dir
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

optlist,alist=getopt.getopt(sys.argv[1:],'hc:m:o:')
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-c':count_path = opt[1]
    elif opt[0] == '-m':meta_path = opt[1]
    elif opt[0] == '-o':tool_output_dir = opt[1]


min_cell = 10


if __name__ == '__main__':
    
    if tool_output_dir[-1] == '/':
        tool_output_dir = tool_output_dir[:-1]
    prepare_output_path = '{}/CytoTalk/input/'.format(tool_output_dir)
    if not os.path.exists(prepare_output_path):
        os.mkdir(prepare_output_path)
 
    count_df = pd.read_csv(count_path, sep='\t', index_col = 0)

    meta_df = pd.read_csv(meta_path, sep='\t', header=None)
    meta_df.iloc[:,1] = [ct.replace('/','_').replace(' ','_') for ct in meta_df.iloc[:,1]]

    ct_barcode_dic = {}
    for i in range(meta_df.shape[0]):
        if meta_df.iloc[i,1] not in ct_barcode_dic:
            ct_barcode_dic.update({meta_df.iloc[i,1]:[meta_df.iloc[i,0]]})
        else:
            ct_barcode_dic[meta_df.iloc[i,1]].append(meta_df.iloc[i,0])

    for ct in list(ct_barcode_dic.keys()):
        if len(ct_barcode_dic[ct]) < min_cell:
            del(ct_barcode_dic[ct])

    print(ct_barcode_dic.keys())
    for ct in ct_barcode_dic.keys():
        ct_list = ct_barcode_dic[ct]
        select_col = [col in ct_list for col in count_df.columns]
        ct_df = count_df.loc[:,select_col]
        ct_df.to_csv('{}scRNAseq_{}.csv'.format(prepare_output_path, ct))
        
    os.system('cp /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/CytoTalk/PAAD_st/Input/deter_W4.R {}'.format(prepare_output_path))
    os.system('cp /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/CytoTalk/PAAD_st/Input/*.txt {}'.format(prepare_output_path))
    os.system('cp /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/CytoTalk/PAAD_st/Input/gen_PCSF.py {}'.format(prepare_output_path))
    
    ### generate the cytotalk.input
    cytotalk_submit_input = '{}/CytoTalk/script/cytotalk.input'.format(tool_output_dir)
    
    # ct_combine = [it for it in itertools.product(list(ct_barcode_dic.keys()),list(ct_barcode_dic.keys())) if it[0] != it[1]]
    ct_combine = list(itertools.combinations(list(ct_barcode_dic.keys()), 2))
    
    input_df = pd.DataFrame(ct_combine)
    input_df['input'] = prepare_output_path
    input_df['output'] = '{}/CytoTalk/output/'.format(tool_output_dir)
    nohup_dir = '{}/CytoTalk/script/nohup'.format(tool_output_dir)
    if not os.path.exists(nohup_dir):
        os.mkdir(nohup_dir)
    input_df['nohup'] = input_df.apply(lambda x: '{}/{}_{}_nohup.out'.format(nohup_dir,x[0],x[1]), axis=1)
    input_df.to_csv(cytotalk_submit_input, sep='\t', header=None, index=None)
    
    ### generate submit code
    submit_code = '''cat {} | while read ct_a ct_b input_dir output_dir nohup_path
do
while true
do
CURRENT=$(ps -ef | grep liu | grep -E 'time -v Rscript' | grep -E 'cytotalk' | wc -l)
if [ ${{CURRENT}} -lt 3 ];then
screen -x cytotalk -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_cytotalk.R -a $ct_a -b $ct_b -i $input_dir -o $output_dir > $nohup_path 2>&1 &\\n"
sleep 20s
break
fi
done
done'''.format(cytotalk_submit_input)

    with open('{}/CytoTalk/script/run_cytotalk.sh'.format(tool_output_dir),'w') as f:
        f.write(submit_code)
        
    os.system('chmod 777 {}/CytoTalk/script/run_cytotalk.sh'.format(tool_output_dir))

    prepare_generate_code = 'python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_cytotalk_regenerate.py -i {}/CytoTalk/script/cytotalk.input -o {}/CytoTalk/output -c '.format(tool_output_dir,tool_output_dir)
    
    with open('{}/CytoTalk/script/submit_prepare_regenerate.sh'.format(tool_output_dir),'w') as f:
        f.write(prepare_generate_code)

#     regenerate_result_code = '''cat {} | while read ct_a ct_b input_dir output_dir nohup_path
# do
# while true
# do
# CURRENT=$(ps -ef | grep liu | grep run_cytotalk | wc -l)
# if [ ${{CURRENT}} -lt 6 ];then
# screen -x cytotalk -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/regenerate_cytotalk_result.R -a $ct_a -b $ct_b -i $input_dir -o $output_dir >> $nohup_path 2>&1 &\\n"
# sleep 20s
# break
# fi
# done
# done'''.format(cytotalk_submit_input)

#     with open('{}/CytoTalk/script/regenerate_result_cytotalk.sh'.format(tool_output_dir),'w') as f:
#         f.write(submit_code)
        
#     os.system('chmod 777 {}/CytoTalk/script/regenerate_result_cytotalk.sh'.format(tool_output_dir))
        
        
    submit_regenerate_code = '''cat {}/CytoTalk/script/cytotalk.part.script | while read script_path
do
while true
do
CURRENT=$(ps -ef | grep liu | grep 'Rscript regenerate' | wc -l)
if [ ${{CURRENT}} -lt 3 ];then
screen -x cytotalk -p 0 -X stuff "cd {}/CytoTalk/script/reg_part_script && nohup time -v Rscript $script_path >> {}/CytoTalk/script/regenerate_nohup.out 2>&1 &\n"
sleep 20s
break
fi
done
done'''.format(tool_output_dir,tool_output_dir,tool_output_dir)
    
    with open('{}/CytoTalk/script/submit_regenerate.sh'.format(tool_output_dir),'w') as f:
        f.write(submit_regenerate_code)
        
    os.system('chmod 777 {}/CytoTalk/script/submit_regenerate.sh'.format(tool_output_dir))
        







