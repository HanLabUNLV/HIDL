
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Align 
from Bio.Align import substitution_matrices
import json


# read in protein IDs
dict_ensembl2pdb = {}
with open("ensemblID.linkers.pdb.ss.txt", "r") as file:
  for line in file:
    key, value = line.strip().split(" ")
    try:
      dict_ensembl2pdb[key].add(value)
    except KeyError:
      dict_ensembl2pdb[key] = set()
      dict_ensembl2pdb[key].add(value)
dict_tree = {}
dict_linker = {}
with open("protein2linkernss.txt", "r") as file:
  for line in file:
    prot, genetree, linker = line.strip().split(" ")
    try:
      dict_tree[genetree].add(linker)
    except KeyError:
      dict_tree[genetree] = set()
      dict_tree[genetree].add(linker)
    try:
      dict_linker[linker].add(prot)
    except KeyError:
      dict_linker[linker] = set()
      dict_linker[linker].add(prot)



#alignments = pairwise2.align.localds("LFEFKGEDLTEEEDGGIIRRIQTRGEGY", "GAPLPMEGVDISPKQDEGVLKVIKREGTGTEMPMIGDRVFVHYTGWLLDGTKFDSSLDRKDKFSFDLGKGEVIKAWDIAIATMKVGEVCHITCKPEYAYGSAGSPPKIPPNATLVFEVELFEFKGEDLTEEEDGGIIRRIQTRGEGYAKPNEGAIVEVALEGYYKDKLFDQRELRFEIGEGENLDLPYGLERAIQRMEKGEHSIVYLKPSYAFGSVGKEKFQIPPNAELKYELHLKSFEKAKESWE", matlist.blosum62, -10, -0.5)
#print(alignments[0])

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


linkerdir = "LinkerFastasKCOfiltered/"
tmpdir="/tmp/"
ss_fasta = "ss_seqs.fasta"
pdb_fasta =  SeqIO.to_dict(SeqIO.parse(ss_fasta, "fasta"))

linker_prot_best_struct_alns = {} 
for tree in dict_tree:
  #print (tree)
  for linker in dict_tree[tree]: 
    #print ("\t"+linker)
    linker_prot_best_struct_alns[linker] = {}
    linker_fasta =  SeqIO.to_dict(SeqIO.parse(linkerdir+linker+".fasta.fas", "fasta"))
    for prot in dict_linker[linker]:
      #print ("\t\t"+prot) 
      linker_prot_best_struct_alns[linker][prot] = {}
      linker_seq = linker_fasta[prot].seq.ungap("-")
      structs = dict_ensembl2pdb[prot]
      #ss_seqs = []
      struct_bestscore = 0
      struct_bestpdb = ''
      struct_bestaln = ''
      for pdb in structs:
        #print ("\t\t\t"+pdb) 
        # collect pdb sequence into file
        #ss_seqs.append(pdb_fasta[pdb+":sequence"]) 
        pdb_seq = pdb_fasta[pdb+":sequence"].seq.ungap("-")
        #pdb_seq = pdb_fasta[pdb+":sequence"].seq
        # align linker to pdb
        #alignments = aligner.align(tmpdir+prot+linker+".fas", tmpdir+prot+linker+"ss.fas")
        alignments = pairwise2.align.localds(linker_seq, pdb_seq, matlist.blosum62, -10, -0.5)
        #alignments = pairwise2.align.localms(linker_seq, pdb_seq, 2, -1, -10, -0.5)
        #print(len(alignments))
        if alignments: 
          sortedalns = sorted(alignments, key=lambda x: x.score, reverse=True)
          if (struct_bestscore < sortedalns[0].score): 
            struct_bestpdb = pdb 
            struct_bestscore = sortedalns[0].score
            struct_bestaln = sortedalns[0]
      if (struct_bestscore > 0):
        linker_prot_best_struct_alns[linker][prot][struct_bestpdb] = struct_bestaln
        print(linker)
        print(prot)
        print(struct_bestpdb)
        print(pairwise2.format_alignment(*struct_bestaln))
        exit()
      #print (struct_bestpdb)
      #print (struct_bestaln)
  #print (linker_prot_best_struct_alns)

print (linker_prot_best_struct_alns)
with open('result.json', 'w') as fp:
    json.dump(linker_prot_best_struct_alns, fp)

