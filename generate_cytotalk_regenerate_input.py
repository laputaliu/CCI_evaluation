import pickle as pkl

with open('/fs/home/liuzhaoyang/project/cci_evaluation/mouse_brain/data/pkl/ct_distype_brain.pkl','rb') as f:
    brain_ct_distype = pkl.load(f)
    
brain_ct_list = list(brain_ct_distype_sr.index)
brain_output_list = ['{}\t{}'.format(brain_ct_list[n], 'near' if ct_f == 0 else 'far') for n, ct_f in enumerate(brain_ct_distype_sr) if ct_f != 2]

brain_kept_ct_list = []
with open('/fs/home/liuzhaoyang/project/cci_evaluation/mouse_brain/data/brain_nf_ct_distype.tsv','r') as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        brain_kept_ct_list.append(line[0])

brain_kept_ct_list += ['{}|{}'.format(ct.split('|')[1], ct.split('|')[0]) for ct in brain_kept_ct_list]
brain_kept_ct_list = [ct.replace('|','_') for ct in brain_kept_ct_list]

raw_input_df = pd.read_csv('/fs/home/liuzhaoyang/project/cci_evaluation/mouse_brain/tools/CytoTalk/script/cytotalk.input', sep='\t', header=None)


input_ct_list = ['{}_{}'.format(raw_input_df.iloc[n,0], raw_input_df.iloc[n,1]) for n in range(raw_input_df.shape[0])]
select_index = [ct in brain_kept_ct_list for ct in input_ct_list]
input_df = raw_input_df.loc[select_index,:]
input_df.to_csv('/fs/home/liuzhaoyang/project/cci_evaluation/mouse_brain/tools/CytoTalk/script/cytotalk.input.kept',sep='\t',header=None, index=None)

