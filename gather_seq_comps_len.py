from os import listdir as ls

out = open('linker_seq_comps_lens.tsv', 'w')
out.write('\t'.join(['AA','Linker','length','Proportion\n']))
def readFasta(fn):
    fasta = {}
    for line in fn:
        if line[0]=='>':
            fasta[line[1:-1]] = ''
            key = line[1:-1]
        else:
            fasta[key] += line.strip()
    return fasta

aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
linker_comp = {'total':0}
domain_comp = {'total':0}
for x in aas:
    linker_comp[x] = 0
    domain_comp[x] = 0


for fn in ls('LinkerFastasNewest'):
    if 'fasta.fas' in fn:
            fasta = readFasta(open('LinkerFastasNewest/'+fn,'r'))
	    ID = fn.split('.')[0]
            for taxa in fasta:
                for aa in aas:
                    length = len(fasta[taxa].strip().replace('-',''))
                    if length>=1:
                        proportion = float(fasta[taxa].count(aa))/length
                        out.write('\t'.join([aa,ID,str(length),str(proportion)])+'\n')
                    


