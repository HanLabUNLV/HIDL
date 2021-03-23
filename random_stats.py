from os import listdir as ls
from statistics import pstdev as sd

#read in fasta file, gonna need to do a lot of this
def readFasta(fn):
    fasta = {}
    for line in fn:
        if line[0]=='>':
            fasta[line[1:-1]] = ''
            key = line[1:-1]
        else:
            fasta[key] += line.strip()
    return fasta
    
#iterate over fastas in a directory
d = 'LinkerFastasNew'
gl_means = []
gl_sds = []
gl_gap_means = []
gl_gap_sds = []
for fn in ls(d):
    fasta = readFasta(open(d+'/'+fn,'r'))
    lengths = []
    gaps = []
    for entry in fasta:
        lengths.append(len(fasta[entry]))
        gaps.append(fasta[entry].count('-'))
    gl_means.append(str(sum(lengths)/len(lengths))+'\n')
    gl_sds.append(str(sd(lengths))+'\n')
    gl_gap_means.append(str(sum(gaps)/len(lengths))+'\n')
    gl_gap_sds.append(str(sd(gaps))+'\n')

out_gl_means = open('Stats/'+d+'_means.txt','w')
out_gl_means.writelines(gl_means)
out_gl_means = open('Stats/'+d+'_sds.txt','w')
out_gl_means.writelines(gl_sds)
out_gl_means = open('Stats/'+d+'_gap_means.txt','w')
out_gl_means.writelines(gl_gap_means)
out_gl_means = open('Stats/'+d+'_gap_sds.txt','w')
out_gl_means.writelines(gl_gap_sds)

