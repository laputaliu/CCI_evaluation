# with open('/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools/CytoTalk/script/ct_kept_list.txt','w') as f:
#     f.write('\n'.join(reoutput_list))


import sys
import os
import getopt

usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -o output dir
    -a available output dir
    -c kept cell tyep list file
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

optlist,alist=getopt.getopt(sys.argv[1:],'ho:a:c:')
for opt in optlist:
    if opt[0] == '-h':sys.exit(usage)
    elif opt[0] == '-o':output_dir = opt[1]
    elif opt[0] == '-a':avoutput_path = opt[1]
    elif opt[0] == '-c':ct_kept_path = opt[1]



if __name__ == '__main__':

    with open(ct_kept_path,'r') as f:
        reoutput_list = [line.strip() for line in f.readlines()]
#         f.write('\n'.join(reoutput_list))

#     output_dir = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools/CytoTalk/output/'
#     avoutput_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_SCC/tools/CytoTalk/output_available'
    if not os.path.exists(avoutput_path):
        os.mkdir(avoutput_path)
    for re_dir in reoutput_list:
        tmp_dir_list = [d.strip() for d in os.popen('ls {}/{}'.format(output_dir, re_dir))]
        if 'IllustratePCSF' in tmp_dir_list:
            if 'PCSF_CrosstalkEdge.txt' in [f.strip() for f in os.popen('ls {}/{}/IllustratePCSF'.format(output_dir, re_dir))]:
                if not os.path.exists('{}/{}/IllustratePCSF'.format(avoutput_path,re_dir)):
                    os.system('mkdir -p {}/{}/IllustratePCSF'.format(avoutput_path,re_dir))
                os.system('cp {}/{}/IllustratePCSF/PCSF_CrosstalkEdge.txt {}/{}/IllustratePCSF/PCSF_CrosstalkEdge.txt'.format(output_dir, re_dir, avoutput_path, re_dir))

