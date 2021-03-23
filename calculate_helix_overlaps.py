#This script will take the files below and calculate the proportion that each overlaps the other
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
out_df = pd.DataFrame(columns=['linker','pid','pdb_slice_prop','linker_slice_prop', 'dssp_slice_prop','filterout'])


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
            origseq, pdbseq, score, start, end = ens_to_pdb[linker][pid][pdbid]
            origseq = origseq.replace('-','')
            pdbseq = pdbseq[start:end]

            if pdbid in pdb:
                    newstart = pdb[pdbid].index(pdbseq.replace('-',''))
                    newend = newstart + (end-start)

                    pdbstrip = pdbseq.strip('-')
                    dsspseq = dssp[pdbid][newstart:newend]
                    
                    #gather linker, pid
                    data = [linker, pid]
                    
                    #proportion of total pdb that this aligned slice constitutes
                    data.append(float(len(pdbstrip))/len(pdb[pdbid]))

                    #proportion of pdb slice that covers linker
                    data.append(float(len(pdbseq))/len(origseq))

                    '''
                    if float(len(pdbseq))/len(origseq) >1:
                        print(ens_to_pdb[linker][pid][pdbid])
                    '''

                    #proportion of non-gap chars in dssp to length of original seq
                    data.append(float(len(dsspseq.replace('-','')))/len(origseq))

                    #filter out bool based on coverage
                    data.append( (float(len(pdbseq))/len(origseq)) < 0.9)

                    row = pd.DataFrame([data], columns=['linker','pid','pdb_slice_prop','linker_slice_prop', 'dssp_slice_prop','filterout'])
                    out_df = out_df.append(row, ignore_index=True)

print('printing to file helix_overlaps.csv...')
out_df.to_csv('helix_overlaps.csv', index=False)




