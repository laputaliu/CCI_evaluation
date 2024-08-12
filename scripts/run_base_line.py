import pandas as pd
import numpy as np
import itertools
import os
import sys
import pickle as pkl
import random

import getopt

usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -c sc count matrix path
    -m sc meta path
    -n ntop_ip
    -o tool output dir
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

# set default
ntop_ip = 20

optlist, alist=getopt.getopt(sys.argv[1:],
                            'hc:m:n:o:',
                            [
                                'help=',
                                'sc_norm=',
                                'sc_meta=',
                                'ntop=',
                                'output_dir='
                            ])
for opt, arg in optlist:
    if opt in ['-h', '--help']:
        sys.exit(usage)
    elif opt in ['-c', '--sc_norm']:
        sc_count_path = arg
    elif opt in ['-m', '--sc_meta']:
        sc_meta_path = arg
    elif opt in ['-n', '--ntop']:
        ntop_ip = int(arg)
    elif opt in ['-o', '--output_dir']:
        output_dir = arg


        
def run_base_line(sc_count_path, sc_meta_path, ntop_ip, output_dir):
    
    cc_ip_path = '../cci_database/cc_ip_all_multi_split_only_single_uniq.tsv'

    sc_count_df = pd.read_csv(sc_count_path, sep='\t', index_col=0)
    sc_meta_df = pd.read_csv(sc_meta_path, sep='\t', index_col=0, header=None)
    sc_meta_df.columns = ['celltype']
    cc_ip_df = pd.read_csv(cc_ip_path, sep='\t', header=None)
    cc_ip_df.columns = ['Ligand', 'Receptor']
    cc_ip_df['ip'] = cc_ip_df.apply(lambda x: '{} - {}'.format(x[0], x[1]), axis=1)
    
    with open('./evaluation_result/pkl/tool_res_dic_direct.pkl', 'rb') as f:
        ori_res_dic = pkl.load(f)
        
    print(ori_res_dic['italk'].keys())

    all_cts = list(set(sc_meta_df.loc[:,'celltype']))
    all_genes = list(set(sc_count_df.index))

    # filter avi lig/rec genes in count data
    lig_genes = list(set(cc_ip_df.loc[:,'Ligand']))
    avi_lig_genes = list(set(lig_genes).intersection(set(all_genes)))
    rec_genes = list(set(cc_ip_df.loc[:,'Receptor']))
    avi_rec_genes = list(set(rec_genes).intersection(set(all_genes)))
    lr_genes = list(set(lig_genes + rec_genes))
    avi_lr_genes = list(set(avi_lig_genes + avi_rec_genes))

    # filter avi lr pairs
    select_index = [(cc_ip_df.iloc[i,0] in avi_lig_genes and cc_ip_df.iloc[i,1] in avi_rec_genes) for i in range(cc_ip_df.shape[0])]
    avi_cc_ip_df = cc_ip_df.loc[select_index,:]
    avi_cc_ip_df.loc[:,'ip'] = avi_cc_ip_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis=1)
    avi_ips = list(set(avi_cc_ip_df.loc[:,'ip']))

    sc_count_df = sc_count_df.loc[avi_lr_genes,:]

    ###################
    # generate result #
    ###################
    
    tool_res_dic = {}

    for ctp in itertools.permutations(all_cts, 2):
        ct_a, ct_b = ctp
        ct_a_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_a,:].index)
        ct_b_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_b,:].index)

        ct_a_count_df = sc_count_df.loc[avi_lig_genes, ct_a_bars]
        ct_b_count_df = sc_count_df.loc[avi_rec_genes, ct_b_bars]

        ct_a_avi_count = ct_a_count_df.apply(np.mean, axis=1)
        ct_b_avi_count = ct_b_count_df.apply(np.mean, axis=1)

        tmp_lig_exp = ct_a_avi_count.loc[list(avi_cc_ip_df.loc[:,'Ligand'])]
        tmp_lig_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
        tmp_rec_exp = ct_b_avi_count.loc[list(avi_cc_ip_df.loc[:,'Receptor'])]
        tmp_rec_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
        tmp_lr_exp = tmp_lig_exp.multiply(tmp_rec_exp)
        tmp_lr_exp = tmp_lr_exp.sort_values(ascending=False)

        if tmp_lr_exp.shape[0] > ntop_ip:
            tool_res_dic.update({
                '{}|{}'.format(ct_a,ct_b): list(tmp_lr_exp[:ntop_ip].index)
            })
        else:
            tool_res_dic.update({
                '{}|{}'.format(ct_a,ct_b): list(tmp_lr_exp.index)
            })



    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(os.path.join(output_dir, 'tool_res_dic.pkl'), 'wb') as f:
        pkl.dump(tool_res_dic, f)

        
# def run_base_line(sc_count_path, sc_meta_path, ntop_ip, output_dir):
    
#     cc_ip_path = '../cci_database/cc_ip_all_multi_split_only_single_uniq.tsv'

#     sc_count_df = pd.read_csv(sc_count_path, sep='\t', index_col=0)
#     sc_meta_df = pd.read_csv(sc_meta_path, sep='\t', index_col=0, header=None)
#     sc_meta_df.columns = ['celltype']
#     cc_ip_df = pd.read_csv(cc_ip_path, sep='\t', header=None)
#     cc_ip_df.columns = ['Ligand', 'Receptor']
#     cc_ip_df['ip'] = cc_ip_df.apply(lambda x: '{} - {}'.format(x[0], x[1]), axis=1)
    
#     with open('./evaluation_result/pkl/tool_res_dic_direct.pkl', 'rb') as f:
#         ori_res_dic = pkl.load(f)
        
#     print(ori_res_dic['italk'].keys())

#     all_cts = list(set(sc_meta_df.loc[:,'celltype']))
#     all_genes = list(set(sc_count_df.index))

#     # filter avi lig/rec genes in count data
#     lig_genes = list(set(cc_ip_df.loc[:,'Ligand']))
#     avi_lig_genes = list(set(lig_genes).intersection(set(all_genes)))
#     rec_genes = list(set(cc_ip_df.loc[:,'Receptor']))
#     avi_rec_genes = list(set(rec_genes).intersection(set(all_genes)))
#     lr_genes = list(set(lig_genes + rec_genes))
#     avi_lr_genes = list(set(avi_lig_genes + avi_rec_genes))

#     # filter avi lr pairs
#     select_index = [(cc_ip_df.iloc[i,0] in avi_lig_genes and cc_ip_df.iloc[i,1] in avi_rec_genes) for i in range(cc_ip_df.shape[0])]
#     avi_cc_ip_df = cc_ip_df.loc[select_index,:]
#     avi_cc_ip_df.loc[:,'ip'] = avi_cc_ip_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis=1)
#     avi_ips = list(set(avi_cc_ip_df.loc[:,'ip']))

#     sc_count_df = sc_count_df.loc[avi_lr_genes,:]

#     ###################
#     # generate result #
#     ###################
    
#     tool_res_dic = {}

#     for ctp in itertools.permutations(all_cts, 2):
#         ct_a, ct_b = ctp
#         ct_a_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_a,:].index)
#         ct_b_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_b,:].index)

#         ct_a_count_df = sc_count_df.loc[avi_lig_genes, ct_a_bars]
#         ct_b_count_df = sc_count_df.loc[avi_rec_genes, ct_b_bars]

#         ct_a_avi_count = ct_a_count_df.apply(np.mean, axis=1)
#         ct_b_avi_count = ct_b_count_df.apply(np.mean, axis=1)

#         tmp_lig_exp = ct_a_avi_count.loc[list(avi_cc_ip_df.loc[:,'Ligand'])]
#         tmp_lig_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
#         tmp_rec_exp = ct_b_avi_count.loc[list(avi_cc_ip_df.loc[:,'Receptor'])]
#         tmp_rec_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
#         tmp_lr_exp = tmp_lig_exp.multiply(tmp_rec_exp)
#         tmp_lr_exp = tmp_lr_exp.sort_values(ascending=False)
        

#         if tmp_lr_exp.shape[0] > ntop_ip:
#             if '{}|{}'.format(ct_a,ct_b) in ori_res_dic['italk']:
#                 print('{}|{} in'.format(ct_a,ct_b))
#                 tmp_ips = []
#                 tmp_ips = random.sample(list(tmp_lr_exp[:50].index), k=3)
#                 it_ips = ori_res_dic['italk']['{}|{}'.format(ct_a,ct_b)]
#                 if len(it_ips) >= ntop_ip-3:
#                     tmp_ips += random.sample(it_ips, k=ntop_ip-3)
#                 else:
#                     tmp_ips += it_ips
#                     tmp_ips += random.sample(list(cc_ip_df.loc[:,'ip']), k=ntop_ip-len(tmp_ips))
#                 tmp_ips = list(set(tmp_ips))
#             else:
#                 print('{}|{} out'.format(ct_a,ct_b))
#                 tmp_ips = []
# #                 tmp_ips = random.sample(list(tmp_lr_exp.index), k=ntop_ip)
            
#             tool_res_dic.update({
#                 '{}|{}'.format(ct_a,ct_b): tmp_ips
#             })
#         else:
#             tool_res_dic.update({
#                 '{}|{}'.format(ct_a,ct_b): list(tmp_lr_exp.index)
#             })



#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     with open(os.path.join(output_dir, 'tool_res_dic.pkl'), 'wb') as f:
#         pkl.dump(tool_res_dic, f)

        

if __name__ == '__main__':
    
    
    run_base_line(sc_count_path, sc_meta_path, ntop_ip, output_dir)
    


    
    
    
    
    
# cc_ip_path = '/fs/home/liuzhaoyang/data/cc_ip/cc_ip_all_multi_split_only_single_uniq.tsv'


# print(sc_count_path)
# print(sc_meta_path)


# sc_count_df = pd.read_csv(sc_count_path, sep='\t', index_col=0)
# sc_meta_df = pd.read_csv(sc_meta_path, sep='\t', index_col=0, header=None)
# sc_meta_df.columns = ['celltype']
# cc_ip_df = pd.read_csv(cc_ip_path, sep='\t', header=None)
# cc_ip_df.columns = ['Ligand', 'Receptor']

# print(sc_meta_df)

# all_cts = list(set(sc_meta_df.loc[:,'celltype']))
# all_genes = list(set(sc_count_df.index))

# print(all_cts)
# print(cc_ip_df)
    
# # filter avi lig/rec genes in count data
# lig_genes = list(set(cc_ip_df.loc[:,'Ligand']))
# avi_lig_genes = list(set(lig_genes).intersection(set(all_genes)))
# rec_genes = list(set(cc_ip_df.loc[:,'Receptor']))
# avi_rec_genes = list(set(rec_genes).intersection(set(all_genes)))
# lr_genes = list(set(lig_genes + rec_genes))
# avi_lr_genes = list(set(avi_lig_genes + avi_rec_genes))

# print(avi_lig_genes)
# print(avi_rec_genes)

# # filter avi lr pairs
# select_index = [(cc_ip_df.iloc[i,0] in avi_lig_genes and cc_ip_df.iloc[i,1] in avi_rec_genes) for i in range(cc_ip_df.shape[0])]
# avi_cc_ip_df = cc_ip_df.loc[select_index,:]
# avi_cc_ip_df.loc[:,'ip'] = avi_cc_ip_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis=1)
# avi_ips = list(set(avi_cc_ip_df.loc[:,'ip']))

# print(avi_cc_ip_df)

# sc_count_df = sc_count_df.loc[avi_lr_genes,:]



# tool_res_dic = {}

# for ctp in itertools.permutations(all_cts, 2):
#     ct_a, ct_b = ctp
#     print(ct_a, ct_b, ctp)
#     ct_a_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_a,:].index)
#     ct_b_bars = (sc_meta_df.loc[sc_meta_df.celltype == ct_b,:].index)
    
#     ct_a_count_df = sc_count_df.loc[avi_lig_genes, ct_a_bars]
#     ct_b_count_df = sc_count_df.loc[avi_rec_genes, ct_b_bars]
    
#     ct_a_avi_count = ct_a_count_df.apply(np.mean, axis=1)
#     ct_b_avi_count = ct_b_count_df.apply(np.mean, axis=1)
    
#     tmp_lig_exp = ct_a_avi_count.loc[list(avi_cc_ip_df.loc[:,'Ligand'])]
#     tmp_lig_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
#     tmp_rec_exp = ct_b_avi_count.loc[list(avi_cc_ip_df.loc[:,'Receptor'])]
#     tmp_rec_exp.index = list(avi_cc_ip_df.loc[:,'ip'])
#     tmp_lr_exp = tmp_lig_exp.multiply(tmp_rec_exp)
    
#     if '{}|{}'.format(ct_a,ct_b) == 'mDCs|Endothelial':
#         print(tmp_lig_exp)
#         print(tmp_rec_exp)
#         print(tmp_lr_exp)
        
#     tmp_lr_exp = tmp_lr_exp.sort_values(ascending=False)
    
#     if '{}|{}'.format(ct_a,ct_b) == 'mDCs|Endothelial':
#         print(tmp_lr_exp)
    
#     if tmp_lr_exp.shape[0] > ntop_ip:
#         tool_res_dic.update({
#             '{}|{}'.format(ct_a,ct_b): list(tmp_lr_exp[:ntop_ip].index)
#         })
#     else:
#         tool_res_dic.update({
#             '{}|{}'.format(ct_a,ct_b): list(tmp_lr_exp.index)
#         })
        

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# with open(os.path.join(output_dir, 'tool_res_dic.pkl'), 'wb') as f:
#     pkl.dump(tool_res_dic, f)
    
