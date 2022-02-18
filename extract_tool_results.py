import numpy as np
import pandas as pd
import os
import pickle as pkl
import re


def extract_ip_cc(pval_path, ct_list, name_trans_path, ct_name_correct=None, cc_ip_path = None, ntop=None):
    
    res_dic = {}
    ## key: ct;
    ## value: [ip_pair(str)]
    
    ip_pval_df = pd.read_csv(pval_path, sep = '\t')
    
    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])
    
    if ct_name_correct != None:        
        tmp_col = list(ip_pval_df.columns)
        for raw_name in ct_name_correct:
            tmp_col = [ct.replace(raw_name, ct_name_correct[raw_name]) for ct in tmp_col]
        ip_pval_df.columns = tmp_col
        
    ip_name_trans_dic = {}
    with open(name_trans_path,'r') as f:
        for line in f.readlines()[1:]:
            ip_name_trans_dic.update({line.strip().split('\t')[0]:line.strip().split('\t')[1]})


        
    for ct_type in ct_list:
        if ntop != None:
            tmp_pval_df = pd.DataFrame({'ip':[],'pval':[]})
        tmp_ip_list = []
        for ct in ['{}|{}'.format(ct_type.split('|')[1],ct_type.split('|')[0]), ct_type]:
            if ct not in list(ip_pval_df.columns):
                continue
            ct_pval_df = ip_pval_df.loc[:,ct]
            ct_pval_df = ct_pval_df.loc[ct_pval_df < 0.05]
            ### ct_pval_df to store the significant ips.
            
            if ntop != None:
                tmp_pval_df = tmp_pval_df.append(pd.DataFrame({'ip':list(ct_pval_df.index),'pval':list(ct_pval_df)}))
                
                
            sig_ips = list(ct_pval_df.index)
                
            sig_ips = [ip_name_trans_dic[ip] for ip in sig_ips]            
            tmp_ip_list += sig_ips
            
        if ntop != None:
            tmp_pval_df = tmp_pval_df.drop_duplicates()
            tmp_pval_df = tmp_pval_df.sort_values(by='pval',ascending=True)
            tmp_ip_list = [ip_name_trans_dic[ip] for ip in list(tmp_pval_df['ip'])]
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))       
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            
        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_cpdb(pval_path, ct_list, ip_trans_path, dec_path, ct_name_correct = None, cc_ip_path = None, ntop=None):
    
    res_dic = {}
    ## key: ct;
    ## value: [ip_pair(str)]

    ip_pval_df = pd.read_csv(pval_path,sep= '\t',index_col = 1)
    
    ip_name_trans_dic = {}
    with open(ip_trans_path,'r') as f:
        for line in f.readlines()[1:]:
            ip_name_trans_dic.update({line.strip().split('\t')[0]:line.strip().split('\t')[1]})

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])

    if ct_name_correct != None:
        tmp_col = list(ip_pval_df.columns)
        for raw_name in ct_name_correct:
            tmp_col = [ct.replace(raw_name, ct_name_correct[raw_name]) for ct in tmp_col]
        ip_pval_df.columns = tmp_col
        
    for ct_type in ct_list:
        tmp_ip_list = []
        if ntop != None:
            tmp_pval_df = pd.DataFrame({'ip':[],'pval':[]})

        for ct in ['{}|{}'.format(ct_type.split('|')[1],ct_type.split('|')[0]), ct_type]:
            if ct not in list(ip_pval_df.columns):
                continue
        
            ct_pval_df = ip_pval_df.loc[:,ct]
            ct_pval_df = ct_pval_df.loc[ct_pval_df < 0.05]                

            cpdb_name_trans_dic = {}
            with open(dec_path,'r') as f:
                for line in f.readlines()[1:]:
                    cpdb_name_trans_dic.update({line.strip().split('\t')[1]:line.strip().split('\t')[0]})
            new_cpdb_index = []
            for ip in list(ct_pval_df.index):
                tmp_list = []
                for gene in ip.split('_'):
                    if gene in cpdb_name_trans_dic:
                        gene = cpdb_name_trans_dic[gene]
                    tmp_list.append(gene)
                new_cpdb_index.append('_'.join(tmp_list))
            ct_pval_df.index = new_cpdb_index
            
            if ntop != None:
                tmp_pval_df = tmp_pval_df.append(pd.DataFrame({'ip':list(ct_pval_df.index),'pval':list(ct_pval_df)}))
            
            tmp_ip_list += [ip_name_trans_dic[ip] for ip in list(ct_pval_df.index) if ip in ip_name_trans_dic]
#         tmp_ip_list = list(set(tmp_ip_list))
        
        if ntop != None:
            tmp_pval_df = tmp_pval_df.drop_duplicates()
            tmp_pval_df = tmp_pval_df.sort_values(by='pval',ascending=True)
            tmp_ip_list = [ip_name_trans_dic[ip] for ip in list(tmp_pval_df['ip']) if ip in ip_name_trans_dic]
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))       
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

#         if cc_ip_path != None:
#             tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
#         if ntop != None:
#             tmp_ip_list = tmp_ip_list[:ntop]

        res_dic.update({ct_type:tmp_ip_list})
    
    return res_dic


def extract_ip_italk(res_path,ct_list,ct_name_correct=None, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part

    res_dic = {}
    
    res_df = pd.read_csv(res_path,sep='\t')
    ## filter out 0 exp
    res_df = res_df.loc[res_df.cell_from_mean_exprs > 0,:]
    res_df = res_df.loc[res_df.cell_to_mean_exprs > 0, :]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])

    if ct_name_correct != None:
        res_df['cell_from'] = res_df['cell_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['cell_to'] = res_df['cell_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        if ct_a in set(res_df['cell_from']) and ct_b in set(res_df['cell_to']):
            sub_df = res_df.loc[res_df.cell_from == ct_a,:].loc[res_df.cell_to == ct_b,:]
            sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
            tmp_ip_list += list(sub_df['ip'])
        if ct_b in set(res_df['cell_from']) and ct_a in set(res_df['cell_to']):
            sub_df = res_df.loc[res_df.cell_from == ct_b,:].loc[res_df.cell_to == ct_a,:]
            sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
            tmp_ip_list += list(sub_df['ip'])
            
        if ntop != None:
            if ct_a in set(res_df['cell_from']) and ct_b in set(res_df['cell_to']) and ct_b in set(res_df['cell_from']) and ct_a in set(res_df['cell_to']):
                tmp_score_df = res_df.loc[res_df.cell_from == ct_b,:].loc[res_df.cell_to == ct_a,:]
                tmp_score_df = tmp_score_df.append(res_df.loc[res_df.cell_from == ct_a,:].loc[res_df.cell_to == ct_b,:])
                tmp_score_df['ip'] = tmp_score_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
                tmp_score_df['score'] = tmp_score_df['cell_from_mean_exprs'] * tmp_score_df['cell_to_mean_exprs']
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:         
            tmp_ip_list = list(set(tmp_ip_list))   
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_scr(res_dir,ct_list, lr_score=0.4, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    
    res_dic = {}
    res_ct_list = [f.split('.')[0] for f in os.listdir(res_dir)]
    ## the output result file name list
    ## result file name format : Ductal_Endothelial.tsv

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        # ct_a for sender celltype
        # ct_b for recevier celltype
        # ct_a -> ct_b
        
        # read scr result file named by celltype pair            
        if '{}_{}'.format(ct_a, ct_b) in res_ct_list and '{}_{}'.format(ct_b, ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for pair_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):

                        ####### extract predicting score. from tools result #######'
                sub_df = pd.read_csv('{}/{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]),sep='\t')
                sub_df = sub_df.loc[sub_df.LRscore > lr_score,:]

                ## the scores of italk result are always larger than 0
                sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
                
                if ntop != None:
                    if pair_index == 0:
                        tmp_score_df = sub_df
                    else:
                        tmp_score_df = tmp_score_df.append(sub_df)

                tmp_ip_list += list(sub_df['ip'])
                
            if ntop != None:
                tmp_score_df = tmp_score_df.drop_duplicates()
                tmp_score_df = tmp_score_df.sort_values(by='LRscore',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
            else:
                tmp_ip_list = list(set(tmp_ip_list))        
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_cytotalk(res_dir,ct_list, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    
    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})

        # read scr result file named by celltype pair            

        if '{}_{}'.format(ct_a, ct_b) in res_ct_list or '{}_{}'.format(ct_b, ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                ct_pair2 = ct_pair.replace('|','_')
                if ct_pair2 in res_ct_list:
                    file_ct = ct_pair2
                    try: 
                        sub_df = pd.read_csv('{}/{}/IllustratePCSF/PCSF_CrosstalkEdge.txt'.format(res_dir, file_ct),sep='\t')
                        sub_df = sub_df.loc[sub_df.CrosstalkScore > 0,:]
                        sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0].split('(')[0].strip(), x[1].split('(')[0].strip()), axis = 1)

                        if ntop != None:
                            tmp_score_df = tmp_score_df.append(sub_df)
                        tmp_ip_list += list(sub_df['ip'])
                    except:
                        continue
        
        if ntop != None and len(tmp_ip_list) > 0:
            tmp_score_df = tmp_score_df.drop_duplicates()
            tmp_score_df = tmp_score_df.sort_values(by='CrosstalkScore',ascending=False)
            tmp_ip_list = list(tmp_score_df['ip'])
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                    
    return res_dic


def extract_ip_natmi(res_path,ct_list, ct_name_correct = None, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    res_dic = {}
    res_df = pd.read_csv(res_path)
    res_df = res_df.loc[:,['Sending cluster','Ligand symbol','Receptor symbol','Target cluster','Edge average expression derived specificity']]
    res_df.columns = ['ct_from','ligand','receptor','ct_to','weight']
    res_df = res_df.loc[res_df.weight >= 0.1,:]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    if ct_name_correct != None:
        res_df['ct_from'] = res_df['ct_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['ct_to'] = res_df['ct_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)

    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        # ct_a for sender celltype
        # ct_b for recevier celltype
        # ct_a -> ct_b

        # read scr result file named by celltype pair            
        if ct_a in set(res_df.ct_from) and ct_b in set(res_df.ct_to):
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                sub_df = res_df.loc[res_df.ct_from == ct_pair.split('|')[0],:].loc[res_df.ct_to == ct_pair.split('|')[1],:]
                if sub_df.shape[0] == 0:
                    continue

                sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[1], x[2]), axis = 1)
                
                if ntop != None:
                    if ct_index == 0:
                        tmp_score_df = sub_df
                    else:
                        tmp_score_df = tmp_score_df.append(sub_df)
                
                tmp_ip_list += list(sub_df['ip'])
        
        if ntop != None and ct_a in set(res_df.ct_from) and ct_b in set(res_df.ct_to):
            tmp_score_df = tmp_score_df.drop_duplicates()
            tmp_score_df = tmp_score_df.sort_values(by='weight',ascending=False)
            tmp_ip_list = list(tmp_score_df['ip'])
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
        
        res_dic.update({ct_type:tmp_ip_list})

                
    return res_dic


def trans_symbol(ip):
    ## be used in extract_ip_icellnet
    if ip.count('+') > 0:
        tmp_ip = ip.split(' - ')
        ip = ' - '.join([tmp_ip[0], '({})'.format('+'.join(tmp_ip[-1].split(' + ')))])
    return ip

def extract_ip_icellnet(res_dir,ct_list, lr_score=20, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    res_dic = {} 
    res_ct_list = [f.split('.')[0][6:] for f in os.listdir(res_dir)]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})
            # split df by celltype
        if '{}_{}'.format(ct_a,ct_b) in res_ct_list or '{}_{}'.format(ct_b,ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                if ct_pair in ct_distype_sr.index:
                    sub_df = pd.read_csv('{}/score_{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t', header=None)
                    sub_df.columns = ['ip','score']
                    sub_df.dropna(how='any', inplace=True)
                    sub_df = sub_df.loc[sub_df.score > lr_score,:]

                    sub_df['ip'] = sub_df['ip'].apply(lambda x: x.replace(' / ',' - '))
                    sub_df['ip'] = sub_df['ip'].apply(trans_symbol)
                    
                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(sub_df)

                    tmp_ip_list += list(sub_df['ip'])
        
        if ntop != None and ('{}_{}'.format(ct_a,ct_b) in res_ct_list or '{}_{}'.format(ct_b,ct_a) in res_ct_list):
            tmp_score_df = tmp_score_df.drop_duplicates()
            if tmp_score_df.shape[0] == 0:
                tmp_ip_list = []
            else:
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})

    return res_dic
  
def extract_ip_icellnet_refine(res_dir,ct_list, cc_ip_path = None, lr_score_cutoff = 20,ntop=None,pair_top=100):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    res_dic = {} 
    res_ct_file_list = [f for f in os.listdir(res_dir) if f.startswith('score') and f.count('pdf')==0]
    res_ct_list = [f.replace('score_','').replace('_out_all.tsv','') for f in res_ct_file_list]
    

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for fn, ct_type in enumerate(res_ct_list):
        tmp_ip_list = []
        ct_a = ct_type
        
        pair_df = pd.read_csv('{}/score_{}_out_all.tsv'.format(res_dir,ct_a),sep='\t')
        pair_df = pair_df.iloc[:pair_top,:]
        
        ct_b_list = list(pair_df.columns)
        
        for ct_b in ct_b_list:
                
            if ntop != None:
                tmp_score_df = pd.DataFrame({})
                # split df by celltype
            if '{}|{}'.format(ct_a,ct_b) in ct_list:
                ### confirm the result matrix has the specific celltype interaction type
                ct_pair = '{}|{}'.format(ct_a,ct_b)
                if ct_pair in ct_list:
                    sub_df = pair_df.loc[:,ct_b]
                    sub_df.fillna(0,inplace=True)
                    sub_df = sub_df.loc[sub_df>lr_score_cutoff]
                    sub_ip_list = list(sub_df.index)
                    sub_ip_list = [x.replace(' / ',' - ') for x in sub_ip_list]
                    sub_ip_list = [trans_symbol(x) for x in sub_ip_list]
#                         sub_df = pd.read_csv('{}/score_{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t', header=None)
#                     sub_df.columns = ['ip','score']
#                     sub_df.dropna(how='any', inplace=True)
#                     sub_df = sub_df.loc[sub_df.score > lr_score,:]

#                     sub_df['ip'] = sub_df['ip'].apply(lambda x: x.replace(' / ',' - '))
#                     sub_df['ip'] = sub_df['ip'].apply(trans_symbol)

                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(pd.DataFrame({'ip':sub_ip_list,'score':list(sub_df)}))

                    tmp_ip_list += sub_ip_list
        
            if ntop != None and '{}|{}'.format(ct_a,ct_b) in ct_list:
                tmp_score_df = tmp_score_df.drop_duplicates()
                if tmp_score_df.shape[0] == 0:
                    tmp_ip_list = []
                else:
                    tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                    tmp_ip_list = list(tmp_score_df['ip'])
                    if cc_ip_path != None:
                        tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                    tmp_ip_list = tmp_ip_list[:ntop]
            else:            
                tmp_ip_list = list(set(tmp_ip_list))        
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                    
            for ct_dic_type in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_dic_type in ct_list:
                    if ct_dic_type not in res_dic:
                        res_dic.update({ct_dic_type:tmp_ip_list})
                    else:
                        res_dic[ct_dic_type] += tmp_ip_list

#             res_dic.update({ct_type:tmp_ip_list})
    for ct in ct_list:
        if ct not in res_dic:
            res_dic.update({ct:[]})

    return res_dic



def extract_ip_nichenet(res_dir,ct_list, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_list = []   
            
        if 'output{}_{}'.format(ct_a, ct_b) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_pair in ct_distype_sr.index:
                    sub_df = pd.read_csv('{}/output{}_{}/LR_potential.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]),sep='\t')
                    for l_i in range(sub_df.shape[0]):
                        for c_j in range(sub_df.shape[1]):
                            if sub_df.iloc[l_i,c_j] != 0:
                                tmp_ip_list.append('{} - {}'.format(sub_df.columns[c_j], sub_df.index[l_i]))
                                if ntop != None:
                                    tmp_score_list.append(sub_df.iloc[l_i,c_j])
        
        if ntop != None:
            tmp_score_sr = pd.Series(tmp_score_list, index= tmp_ip_list)
            tmp_score_sr = tmp_score_sr.sort_values(ascending=False)
            tmp_ip_list = list(tmp_score_sr.index)
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})

                
    return res_dic
    

def extract_ip_scmlnet(res_dir,ct_list, cc_ip_path = None, ntop=None):
    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
#         if ntop != None:
#             tmp_score_df = pd.DataFrame({})
        
            # split df by celltype
        for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
            if '{}_{}_RecClus.tsv'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list and '{}_{}_LigClus.tsv'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list and '{}_{}_LR.txt'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
                with open('{}/{}_{}_LR.txt'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), 'r') as f:
                    lr_list = [line.strip() for line in f.readlines()]
                tmp_r_df = pd.read_csv('{}/{}_{}_RecClus.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t')
                tmp_l_df = pd.read_csv('{}/{}_{}_LigClus.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t')
                tmp_ip_list += [lr.replace('_',' - ') for lr in lr_list]
        tmp_ip_list = list(set(tmp_ip_list))
        
        if cc_ip_path != None:
            tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
        if ntop != None:
            tmp_ip_list = tmp_ip_list[:ntop]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic





def extract_ip_connectome(res_path, ct_list, ct_name_correct = None, cc_ip_path = None, ntop=None):
    res_dic = {}
    res_df = pd.read_csv(res_path,sep='\t',header=None)
    res_df.columns = ['edge','score']
    res_df['ct_from'] = res_df['edge'].apply(lambda x: x.split(' - ')[0])
    res_df['ct_to'] = res_df['edge'].apply(lambda x: x.split(' - ')[-1])
    res_df['ip'] = res_df['edge'].apply(lambda x: '{} - {}'.format(x.split(' - ')[1], x.split(' - ')[2]))
    ## filter out 0 exp
    res_df = res_df.loc[res_df.score > 0,:]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    if ct_name_correct != None:
        res_df['ct_from'] = res_df['ct_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['ct_to'] = res_df['ct_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})
        
        if ct_a in set(res_df['ct_from']) and ct_b in set(res_df['ct_to']):
            ### confirm the result matrix has the specific celltype interaction type
            for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_pair in ct_distype_sr.index:

                    sub_df = res_df.loc[res_df.ct_from == ct_pair.split('|')[0],:].loc[res_df.ct_to == ct_pair.split('|')[1],:]
                    if sub_df.shape[0] == 0:
                        continue
                    
                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(sub_df)
                        
                    tmp_ip_list += list(sub_df['ip'])
                    
        
        if ntop != None:
            if tmp_score_df.shape[0] == 0:
                tmp_ip_list = []
            else:
                tmp_score_df = tmp_score_df.drop_duplicates()
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic
                
    
    
 if __name__ == '__main__':
    
    with open('/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/pkl/ct_distype_scc.pkl', 'rb') as f:
        ct_distype_sr = pkl.load(f)

    ct_list = list(ct_distype_sr.index)

    #### tools output path #####

    pval_path_cc = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/cc/output/pval_matrix.tsv'
    name_trans_path_cc = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/cc/output/interaction_name_trans.tsv'

    pval_path_cpdb = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/cpdb/output/pvalues.txt'
    name_trans_path_cpdb = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/cc/output/interaction_name_trans.tsv'
    dec_path_cpdb = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/cpdb/output/deconvoluted.txt'

    res_path_italk = '/Volumes/Backup/data/spatial_compare_data/human_SCC/tools/iTALK/output/LR_result.tsv'

    res_path_scr = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/SCR/output'

    res_dir_cytotalk = '/Users/laputaliu/Data/spatial_compare/CytoTalk/output_cellchatdb'

    res_path_natmi = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/NATMI/output/Edges_LRdb_cellchat.csv'

    res_path_icellnet = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/icellnet/output'

    res_dir_nichenet = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/NicheNet/output'

    res_dir_scmlnet = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/scmlnet/output'

    res_path_connectome = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools_TN_Epithelial/Connectome/output/weight_sc.tsv'


    #### extract tools results ####
    tool_res_dic = {}

    res_cc = extract_ip_cc(pval_path_cc, ct_list, name_trans_path_cc)
    print('cc done')

    res_cpdb = extract_ip_cpdb(pval_path_cpdb, ct_list, name_trans_path_cpdb, dec_path_cpdb)
    print('cpdb done')

    res_italk = extract_ip_italk(res_path_italk,ct_list)
    print('italk done')

    res_scr = extract_ip_scr(res_path_scr,ct_list)
    print('scr done')

    res_cytotalk = extract_ip_cytotalk(res_dir_cytotalk,ct_list)
    print('cytotalk done')

    res_natmi = extract_ip_natmi(res_path_natmi, ct_list)
    print('natmi done')

    res_icellnet = extract_ip_icellnet(res_path_icellnet,ct_list)
    print('icellnet done')

    res_nichenet = extract_ip_nichenet(res_dir_nichenet,ct_list)
    print('nichenet done')

    res_scmlnet = extract_ip_scmlnet(res_dir_scmlnet,ct_list)
    print('scmlnet done')

    res_connectome = extract_ip_connectome(res_path_connectome, ct_list)
    print('connectome done')

    tool_res_dic.update({
        'cc':res_cc,
        'cpdb':res_cpdb,
        'italk':res_italk,
        'scr':res_scr,
        'cytotalk':res_cytotalk,
        'natmi':res_natmi,
        'icellnet':res_icellnet,
        'nichenet':res_nichenet,
        'scmlnet':res_scmlnet,
        'connectome':res_connectome
    })
    
    with open('/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/pkl/tool_res_dic_scc.pkl','wb') as f:
        pkl.dump(tool_res_dic, f)

    
    
    
