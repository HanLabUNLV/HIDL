from os import listdir as ls

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
            for taxa in fasta:
                for aa in aas:
                    linker_comp[aa] += fasta[taxa].count(aa)
                    linker_comp['total'] += fasta[taxa].count(aa)

for fn in ls('DomainFastasNewest'):
    if 'fasta.fas' in fn:
            fasta = readFasta(open('DomainFastasNewest/'+fn,'r'))
            for taxa in fasta:
                for aa in aas:
                    domain_comp[aa] += fasta[taxa].count(aa)
                    domain_comp['total'] += fasta[taxa].count(aa)

print('\t'.join(['AA','Linker','Domain']))
for aa in aas:
    print('\t'.join([aa,str(float(linker_comp[aa])/linker_comp['total']), str(float(domain_comp[aa])/domain_comp['total'])]))
