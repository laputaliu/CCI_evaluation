import numpy as np
import ot
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import random
from functools import partial
from concurrent.futures import ThreadPoolExecutor
import os


def data_prepare(count_path,gene_name_list,min_spot_pct=0.1):
    # filter genes
    count_df = pd.read_csv(count_path, sep = '\t',index_col = 0)
    select_index = [gene in gene_name_list for gene in list(count_df.index)]
    count_df = count_df.loc[select_index,:]
    count_df = count_df.loc[np.sum(count_df > 0,axis=1) >= min_spot_pct*count_df.shape[1],:]

    return count_df


def generate_gene_distribution(gene, count_df):
    '''
    extract gene spatial distribution from ST data
    input:
        - gene : gene (type : str)(e.g. 'gene', 'gene1+gene2')
        - count_df : normalized count matrix (gene*spot) (type : pd.DataFrame)
    output:
        - flag : {0,1} : 0 -> gene in count_df; 1 -> gene not in count_df
        - gene_df : Non-zero gene expression value dataframe with spot name as its index (e.g. 'spot1')
                (type : pd.DataFrame)
                
    spot num < min_spot --> trimed
    '''
    
    spot_list = list(count_df.columns)
    
    # check single/multi genes
    if '+' not in gene:
        # simple
        if gene in list(count_df.index):
            gene_df = count_df.loc[gene,:]
            gene_df = gene_df.loc[gene_df > 0]
            if gene_df.shape[0] == 0:
                return 1, None
            return 0, gene_df
        else:
#             print('{} not in counts df'.format(gene))
            return 1, None
    else:
        # complex (multi-subunits)
        gene = gene.split('+')
        for g in gene:
            if g not in list(count_df.index):
#                 print('{} not in counts df'.format(g))
                return 1, None
        gene_df = count_df.loc[gene[0],:]
        sub_num = len(gene)
        for g in gene[1:]:
            gene_df = gene_df*count_df.loc[g,:]
        gene_df = gene_df.apply(lambda x: np.power(x, 1/sub_num))
        gene_df = gene_df.loc[gene_df > 0]
        if gene_df.shape[0] == 0:
            return 1, None
        return 0, gene_df



def cal_gene_distance(ip, count_df, spot_pos_df, reg = 0.01, iternum = 200, shuffle = False):
    '''
    calculate the spatial distance between two gene spatial distributions (ligand and receptor) using Wasserstein distance.
    input:
        - ip : interaction pair: (ligand - receptor)
        - count_df : normalized count matrix (gene*spot)
        - spot_pos_df : coordinates of each spot
        - reg: regularization term for calculating Wasserstein distance using sinkhorn method
        - iternum: the max iteration number for calculating Wasserstein distance using sinkhorn method
        - shuffle: whether do permutation
    output:
        - gene_distance : the Wasserstein distance between two gene spatial distributions
        
    '''
    gene_a, gene_b = ip.split(' - ')
    gene_a = gene_a.strip().replace('(','').replace(')','')
    gene_b = gene_b.strip().replace('(','').replace(')','')
    
    flag_a, gene_df_a = generate_gene_distribution(gene_a,count_df)
    flag_b, gene_df_b = generate_gene_distribution(gene_b,count_df)
    
    if flag_a + flag_b > 0:
        return None

    if shuffle == False:
        gene_pos_a = np.array(spot_pos_df.loc[[spot in list(gene_df_a.index) for spot in list(spot_pos_df.index)],['row','col']])
        gene_pos_b = np.array(spot_pos_df.loc[[spot in list(gene_df_b.index) for spot in list(spot_pos_df.index)],['row','col']])
    else:
        joined_spot_list = list(spot_pos_df.index)
        
        r_select_spot_a = random.sample(joined_spot_list, gene_df_a.shape[0])
        r_select_spot_b = random.sample(joined_spot_list, gene_df_b.shape[0])
        gene_pos_a = np.array(spot_pos_df.loc[[spot in r_select_spot_a for spot in list(spot_pos_df.index)],['row','col']])
        gene_pos_b = np.array(spot_pos_df.loc[[spot in r_select_spot_b for spot in list(spot_pos_df.index)],['row','col']])

    cost_matrix = ot.dist(gene_pos_a,gene_pos_b,metric='euclidean')
    gene_distance_dir = ot.sinkhorn2(np.array(gene_df_a),np.array(gene_df_b), cost_matrix/cost_matrix.max(), reg = reg, numItermax=iternum)
    
    ### the R-L dir
    cost_matrix = ot.dist(gene_pos_b,gene_pos_a,metric='euclidean')
    gene_distance_rev = ot.sinkhorn2(np.array(gene_df_b),np.array(gene_df_a), cost_matrix/cost_matrix.max(), reg = reg, numItermax=iternum)
    
    gene_distance = (gene_distance_dir+gene_distance_rev)/2

    return gene_distance


def compare_gene_distance(ip,count_df,spot_pos_df,reg = 0.01, iternum = 200, shuffle_num = 500):
    '''
    do permutation for each ligand-receptor pair
    calculate the ratio and p-value
    input:
        - ip : interaction pair: (ligand - receptor)
        - count_df : normalized count matrix (gene*spot)
        - spot_pos_df : coordinates of each spot
        - reg: regularization term for calculating Wasserstein distance using sinkhorn method
        - iternum: the max iteration number for calculating Wasserstein distance using sinkhorn method
        - shuffle_num: permutation times
    output:
        - [ip, dis_ori, dis_shu_avg, dis_r, p_val]

    '''
    dis_ori = cal_gene_distance(ip,count_df,spot_pos_df,reg=reg, iternum=iternum)
    
    if dis_ori == None:
#         print('{} None'.format(ip))
        return [ip, dis_ori, None, None, None]
    dis_shu = []
    for n in range(shuffle_num):
        dis_shu.append(cal_gene_distance(ip, count_df, spot_pos_df, reg=reg, iternum=iternum, shuffle=True))
    dis_shu_avg = np.average(dis_shu)
    
    dis_r = dis_ori/dis_shu_avg
    
    n_count = 0
    if dis_r > 1:
        for d in dis_shu:
            if dis_ori < d:
                n_count += 1
    else:
        for d in dis_shu:
            if dis_ori > d:
                n_count += 1
    p_val = n_count/shuffle_num
    
    return [ip, dis_ori, dis_shu_avg, dis_r, p_val]


def processing(idx_path,
               count_path,
               output_path,
               cc_ip_df,
               cc_gene_name,
               reg=0.01,
                    iternum=500,
                    shuffle_num=500
              ):

    ## spot postion & spatial counts
    cc_ip_list = list(cc_ip_df['interaction_name_2'])
    cc_ip_dic = {}
    for i in range(cc_ip_df.shape[0]):
        cc_ip_dic.update({cc_ip_df.iloc[i,-1]:{'pathway':cc_ip_df.iloc[i,0],'ip_type':cc_ip_df.iloc[i,1]}})

    pos_df = pd.read_csv(idx_path,sep = '\t',index_col = 0)
    count_df = data_prepare(count_path,cc_gene_name)

    pool = ThreadPoolExecutor(4)

    compare_gene_distance_partial = partial(compare_gene_distance, count_df = count_df, spot_pos_df = pos_df, reg=reg, iternum=iternum, shuffle_num=shuffle_num)

    results_iter = pool.map(compare_gene_distance_partial, cc_ip_list)
    res_list_all = [res for res in results_iter if res[-1] != None]
    res_df = pd.DataFrame(res_list_all)
    res_df.columns = ['ip','d_ori','d_shu','d_rat','p_val']
    res_df = res_df.sort_values(by='p_val')
    res_df['pathway'] = res_df.apply(lambda x : cc_ip_dic[x[0]]['pathway'], axis=1)
    res_df['ip_type'] = res_df.apply(lambda x : cc_ip_dic[x[0]]['ip_type'], axis=1)
    
    res_df_long = res_df.loc[res_df.p_val <= 0.1,:].loc[res_df.d_rat > 1]
    res_df_long = res_df_long.sort_values(by='p_val')
    res_df_short = res_df.loc[res_df.p_val <= 0.1,:].loc[res_df.d_rat < 1]
    res_df_short = res_df_short.sort_values(by='p_val')
   
    res_df.to_csv('{}/ip_distance_all.tsv'.format(output_path),sep = '\t', index = None)
    res_df_long.to_csv('{}/ip_distance_long.tsv'.format(output_path),sep = '\t', index = None)
    res_df_short.to_csv('{}/ip_distance_short.tsv'.format(output_path),sep = '\t', index = None)
    
    print('{} done'.format(idx_path.split('/')[-1]))
    
    
if __name__ == '__main__':
  
  ## setting data file path
  st_count_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/data/processed/intestinal_st_sct_norm.tsv'
  st_coord_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/data/processed/intestinal_st_coord.tsv'
  cc_ip_path = '/fs/home/liuzhaoyang/data/cc_ip/cc_ip_all_multi_split.tsv'
  cc_gene_name_path = '/fs/home/liuzhaoyang/project/cci_evaluation/Seurat_mapping/10X_mouse_brain/cc_gene_name_all.txt'
  output_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_heart/data/ip_dis_avgdis_random_all_spot'


  if not os.path.exists(output_path):
      os.mkdir(output_path)

  cc_ip_df = pd.read_csv(cc_ip_path,sep='\t',index_col=0)
  with open(cc_gene_name_path,'r') as f:
      cc_gene_name = [line.strip() for line in f.readlines()]


  processing(st_coord_path,st_count_path,output_path,cc_ip_df,cc_gene_name,shuffle_num=500, reg=0.001)

