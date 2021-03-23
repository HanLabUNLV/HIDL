#USAGE: python remove_duplicate_in_fasta.py filein

from sys import argv

fastafile = open(argv[1],'r')
lines = []
pids = []
for line in fastafile:
    if line in pids:
        break
    elif line[0]=='>':    
        pids.append(line)
        lines.append(line)
    else:
        lines.append(line)
fastafile.close()
fastafile = open(argv[1],'w')
for line in lines:
    fastafile.write(line)
