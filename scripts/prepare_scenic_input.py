import pandas as pd
import getopt
import sys


optlist,alist=getopt.getopt(sys.argv[1:],'hc:o:')
for opt in optlist:
    if opt[0] == '-c':count_path = opt[1]
    elif opt[0] == '-o':output_path = opt[1]


count_df = pd.read_csv(count_path, sep='\t', index_col = 0)

count_df = count_df.T

count_df.to_csv(output_path, sep='\t')
