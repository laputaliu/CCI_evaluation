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
    -t ncores
    -k keep ct_list
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

keep_ct_path = None

optlist,alist=getopt.getopt(sys.argv[1:],'hc:m:o:t:k:')
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-c':count_path = opt[1]
    elif opt[0] == '-m':meta_path = opt[1]
    elif opt[0] == '-o':tool_output_dir = opt[1]
    elif opt[0] == '-t':ncores = opt[1]
    elif opt[0] == '-k':keep_ct_path = opt[1]


min_cell = 10


if __name__ == '__main__':
    
    if tool_output_dir[-1] == '/':
        tool_output_dir = tool_output_dir[:-1]

    meta_df = pd.read_csv(meta_path, sep='\t', header=None)

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
            
    ### generate the scmlnet.input
    scmlnet_submit_input = '{}/scmlnet/script/scmlnet.input'.format(tool_output_dir)
    
    ct_combine = [it for it in itertools.product(list(ct_barcode_dic.keys()),list(ct_barcode_dic.keys())) if it[0] != it[1]]
    
    if keep_ct_path != None:
        with open(keep_ct_path, 'r') as f:
            kept_ct_list = [line.strip() for line in f.readlines()]
        ct_combine_kept = [it for it in ct_combine if '{}_{}'.format(it[0],it[1]) in kept_ct_list or '{}_{}'.format(it[1],it[0]) in kept_ct_list]
        ct_combine = ct_combine_kept
    
    input_df = pd.DataFrame(ct_combine)
    input_df['count'] = count_path
    input_df['meta'] = meta_path
    input_df['ncore'] = ncores
    input_df['output'] = '{}/scmlnet/output/'.format(tool_output_dir)
    nohup_dir = '{}/scmlnet/script/nohup'.format(tool_output_dir)
    if not os.path.exists(nohup_dir):
        os.mkdir(nohup_dir)
    input_df['nohup'] = input_df.apply(lambda x: '{}/{}_{}_nohup.out'.format(nohup_dir,x[0],x[1]), axis=1)
    input_df.to_csv(scmlnet_submit_input, sep='\t', header=None, index=None)
    
    ### generate submit code
    submit_code = '''cat {} | while read ct_a ct_b count_path meta_path ncore output_dir nohup_path
do
while true
do
sleep 60s
CURRENT=$(ps -ef | grep liuzhao | grep scmlnet | wc -l)
if [ ${{CURRENT}} -lt 8 ];then
screen -x scmlnet -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scmlnet_multi.R -a $ct_a -b $ct_b -c $count_path -m $meta_path -n $ncore -o $output_dir > $nohup_path 2>&1 &\\n"
break
fi
done
done'''.format(scmlnet_submit_input)

    with open('{}/scmlnet/script/run_scmlnet.sh'.format(tool_output_dir),'w') as f:
        f.write(submit_code)
        
    os.system('chmod 777 {}/scmlnet/script/run_scmlnet.sh'.format(tool_output_dir))


