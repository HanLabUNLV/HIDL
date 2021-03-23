
from Bio import SeqIO


seqs = [] # Setup an empty list
dssp = [] # Setup an empty list

records = list(SeqIO.parse("ss2.txt", "fasta"))

for i in range(0, len(records)): 
    if i % 2: 
        dssp.append(records[i]) 
    else : 
        seqs.append(records[i]) 

SeqIO.write(seqs, "ss_seqs.fasta", "fasta")
SeqIO.write(dssp, "ss_dssp.fasta", "fasta")

for i in range(0, len(records)): 
    indeces = []
    i = str(record.seq).lower().find(query)
    while i >= 0:
        indeces.append(str(i + 1))
        i =  str(record.seq).lower().find(query, i + 1)

    if len(indeces) > 0:
        print record.id + "\t" + ','.join(indeces)
