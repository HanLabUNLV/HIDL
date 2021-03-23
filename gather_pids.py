#USAGE: python gather_pids.py dir/to/fastas dir/out
from sys import argv
from os import listdir as ls

fastadir = argv[1]
outdir = argv[2]

for fn in ls(fastadir):
    pids = []
    f = open(fastadir + '/' + fn, 'r')
    lines = f.readlines()
    for line in lines:
        if line[0] == '>':
            nextline = lines[lines.index(line)+1]
            if nextline != '\n':
                pids.append(line[1:])        
    out = open(outdir+'/'+fn.split('.')[0]+'_pids.txt','w')
    out.writelines(pids)


