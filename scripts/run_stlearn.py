
import stlearn as st
import numpy as np
import pandas as pd
import sys
import getopt
import itertools


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
	<Required>
	-c st count file
    -p st position (coordinates) file
    -d st deconvolution (STRIDE) result
    -o output dir
    -n random pair number
    -s pval_adj_cutoff'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)        #usage 相当于帮助文件，当无参数调用时输出的帮助文件

# output_path='hw1.3_output.txt'  #设定默认值
# min_length=5
random_pair_num = 1000
pval_adj_cutoff = 0.05
optlist,alist=getopt.getopt(sys.argv[1:],'hc:p:d:o:n:s:') #取出参数
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-c':count_file = opt[1]
    elif opt[0] == '-p':coord_file = opt[1]
    elif opt[0] == '-d':ct_frac_file = opt[1]
    elif opt[0] == '-o':output_dir = opt[1]
    elif opt[0] == '-n':random_pair_num = opt[1]
    elif opt[0] == '-s':pval_adj_cutoff = opt[1]
random_pair_num = int(random_pair_num)
pval_adj_cutoff = float(pval_adj_cutoff)

def write_sig_ip_ct_matrix(adata, target_field, output_dir):
    sig_ip_ct_dic = adata.uns[target_field]
    sig_ips = list(sig_ip_ct_dic.keys())
    sig_celltypes = list(sig_ip_ct_dic[sig_ips[0]].index)
    sig_celltype_pairs = ['|'.join(ctp) for ctp in itertools.permutations(sig_celltypes, 2)]

    sig_ip_ct_count = pd.DataFrame(
        np.zeros((len(sig_ips), len(sig_celltype_pairs))),
        index = sig_ips,
        columns = sig_celltype_pairs
    )

    for tmp_ip in sig_ips:
        for tmp_ctp in sig_celltype_pairs:
            ct_a, ct_b = tmp_ctp.split('|')
            sig_ip_ct_count.loc[tmp_ip,tmp_ctp] = sig_ip_ct_dic[tmp_ip].loc[ct_a,ct_b]

    sig_ip_ct_count.to_csv('{}/{}.tsv'.format(output_dir, target_field), sep='\t')

    
if __name__ == '__main__':
    
    ## actually the stlearn object can be build from 2 pandas.Dataframe
    st_count_df = pd.read_csv(count_file,
                             sep='\t', index_col = 0)
    # the count df need to be row as barcode, column as gene
    st_count_df = st_count_df.T

    # 2 columns are needed (imagecol, imagerow) in st_coord_df
    st_coord_df = pd.read_csv(coord_file,
                             sep='\t', index_col = 0)
    st_coord_df = st_coord_df.loc[:,['row','col']]
    st_coord_df.columns = ['imagerow','imagecol']


    ## create stlearn object using count & spatial information
    adata = st.create_stlearn(count=st_count_df,spatial=st_coord_df,library_id="stlearn")


    ## start stlearn analysis
    adata.var_names_make_unique()
    st.pp.filter_genes(adata, min_cells=3)
    st.pp.normalize_total(adata) # NOTE: no log1p

    # Loading the LR databases available within stlearn (replace them with cellchatDB)
    lrs = []
    with open('/fs/home/liuzhaoyang/data/cc_ip/cc_ip_all_multi_split_only_single.txt','r') as f:
        for line in f.readlines():
            lrs.append(line.strip())
    lrs = list(set(lrs))
    lrs = np.array(lrs)


    # Running the LR analysis #
    st.tl.cci.run(adata, lrs,
                      min_spots = 5, #Filter out any LR pairs with no scores for less than min_spots
                      distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                      n_pairs=random_pair_num, # Number of random pairs to generate; low as example, recommend ~10,000
                      n_cpus=6, # Number of CPUs for parallel. If None, detects & use all available.
                      )

    # p value adjustment
    st.tl.cci.adj_pvals(adata, correct_axis='spot',
                   pval_adj_cutoff=pval_adj_cutoff, adj_method='fdr_bh')

    adata.write('{}/adata_lr_analysis.h5'.format(output_dir))
    
    ## input need dominate (labels) & celltype score (ct_frac)
    ct_frac = pd.read_csv(ct_frac_file,
                         sep = '\t', index_col = 0)

    dominate_ct = [list(ct_frac.columns)[i] for i in ct_frac.apply(np.argmax, axis=1)]

    ## add celltype annotation to data
    adata.obs['cell_type'] = dominate_ct # dominant celltype (only need a list)
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.uns['cell_type'] = ct_frac # cell type fraction

    ## predicting significant CCIs between celltypes
    # Running the counting of co-occurence of cell types and LR expression hotspots #
    st.tl.cci.run_cci(adata, 'cell_type', # Spot cell information either in data.obs or data.uns
                      min_spots=3, # Minimum number of spots for LR to be tested.
                      spot_mixtures=True, # If True will use the label transfer scores,
                                          # so spots can have multiple cell types if score>cell_prop_cutoff
                      cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                      sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                      n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                     )
    
    adata.write('{}/adata_run_cci.h5'.format(output_dir))
    
    
    ## output ip_ct_interaction matrix
    
    for t_f in ['per_lr_cci_cell_type','per_lr_cci_pvals_cell_type','per_lr_cci_raw_cell_type']:
        write_sig_ip_ct_matrix(adata, t_f, output_dir)



