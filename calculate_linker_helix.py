#the goal of this script is to take each linker region, find its mapping onto pdb chain data, then create and use a map between pdb chain data to secondary structure through dssp assignment. The output will be a dataframe with columns linker, PID, and each letter that corresponds to the table quoted below translating dssp output. The value of each of these columns will be the proportion of that linker which is represented by that letter. Of particular interest is the I helix, which is a single amino acid different from the H helix, and are highly detrimental insertions, so their presence may denote functional disorder, especially if paired with high conservation. Don't forget to look for conserved disorder!
'''
    H = alpha-helix
    B = residue in isolated beta-bridge
    E = extended strand, participates in beta ladder
    G = 3-helix (310 helix)
    I = 5 helix (pi-helix)
    T = hydrogen bonded turn
    S = bend 
'''
import pandas as pd
out_df = pd.DataFrame(columns=['linker','pid','H','B','E','G','I','T','S','-','length'])

dssp_chars = ['H','B','E','G','I','T','S','-']

#import multi-line fasta and strip the :sequence or :secstr from the id
def importFasta(fn):
    pid = ''
    seq = ''
    fasta = {}
    for line in fn:
        if '>' in line:
            if seq!='':
                fasta[pid]=seq
                seq=''
            pid = ':'.join(line[1:].split(':')[:-1])
        else:
            seq += line.strip()
    fasta[pid] = seq
    return fasta

print('Opening pdb file...')
#create the pdb chain residue translator, a dictionary with a key that is the name of the sequence in ss_seqs.fasta and a value of the following sequence
fpdb = open('ss_seqs.fasta','r')
print('importing fasta to dict')
pdb = importFasta(fpdb)
fpdb.close()

print('Opening dssp file...')
#create the dssp assignment translator, a dictionary with a key that is the name of the sequence in the ss_dssp.fasta file and a value of the following sequence
fdssp = open('ss_dssp.fasta','r')
dssp = importFasta(fdssp)
fdssp.close()

print('Opening ens_to_pdb file...')
#iterate through linkers in ens_to_pdb. Find first the pdb identifier, then also the snippet of the pdb chain residues that actually applies to the linker. Translate the coordinates of the original linker to match these values so that we can then grab only the relevant portions of the dssp assignment. From these, calculate proportions of each residue in the quoted table above (also '-') and save them in a dataframe
fens = open('result.json','r')
ens_to_pdb = eval(fens.readline())
fens.close()

print('Finding overlap...')
for linker in ens_to_pdb:
    for pid in ens_to_pdb[linker]:
        if len(ens_to_pdb[linker][pid])>1:
            print(ens_to_pdb[linker][pid])
            exit('1:1 pid to pdb not found')
        for pdbid in ens_to_pdb[linker][pid]:
            orig_seq, pdbseq, score, start, end = ens_to_pdb[linker][pid][pdbid]
            origseq = orig_seq[start:end]
            pdbseq = pdbseq[start:end]

            if pdbid in pdb:
                    newstart = pdb[pdbid].index(pdbseq.replace('-',''))
                    newend = newstart + (end-start)

                    dsspseq = dssp[pdbid][newstart:newend]
                    
                    data = [linker, pid]
                    for ch in dssp_chars:
                        data.append(float(dsspseq.count(ch))/len(dsspseq))
                    data.append(len(orig_seq.replace('-','')))
                    row = pd.DataFrame([data], columns=['linker','pid','H','B','E','G','I','T','S','-','length'])
                    out_df = out_df.append(row, ignore_index=True)

print('printing to file...')
out_df.to_csv('helix_proportion_linkers.csv', index=False)




