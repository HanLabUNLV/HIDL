import os
from sys import exit, argv
import re

residue = argv[1]

brown_alphas = []
ou_alphas = []
white_alphas = []
direct = '/data1/home/dbarth/LinkerProject/'+residue+'_residues/'
#fn = '/home/dbarth/Work/GeneTrees/LinkerFastas/model_fit_ENSGT00900000140842_3_Linker_0_1.txt'
for fn in os.listdir(direct):
    if fn[0:5]=='model':
        line_offset=[]
        offset = 0
        aiclist = []
        with open(direct+fn, 'r') as f:
            index = 0
            lines = f.readlines()
            for line in lines:
                if line=='$aic.c\n':
                    aiclist.append(float(lines[index+1].split(' ')[1][:-1]))
                if line=='$alpha\n':
                    alpha = lines[index+2].split(' ')[1][:-1]
                index +=1
            #print alphai
	    try:
		with open(direct+'wh_'+fn,'r') as f:
		    index = 0
		    lines = f.readlines()
		    for line in lines:
			if line=='$aicc\n':
			    aiclist.append(float(lines[index+1].split(' ')[1][:-1]))
			index+=1
	    except:
		pass
	    if len(aiclist)>1:
                if min(aiclist)==aiclist[0]:
                    brown_alphas.append(alpha)
                elif min(aiclist)==aiclist[1]:
                    ou_alphas.append(alpha)
		elif min(aiclist)==aiclist[2]:
		    white_alphas.append(alpha)
with open('/data1/home/dbarth/LinkerProject/'+residue+'_brownian_alphas.txt','w') as f:
    f.write('\n'.join(brown_alphas))
with open('/data1/home/dbarth/LinkerProject/'+residue+'_ou_alphas.txt','w') as f:
    f.write('\n'.join(ou_alphas))
with open('/data1/home/dbarth/LinkerProject/'+residue+'_wh_alphas.txt','w') as f:
    f.write('\n'.join(white_alphas))



