# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:36:00 2018

@author: dylanbarth
"""

from Bio import Phylo
import os

direct = "/home/dbarth/Work/GeneTrees/genePhylos/" #directory with only .xml files
new_files = "/home/dbarth/Work/GeneTrees/geneFastas/" #repository where .fasta goes

for filename in os.listdir(direct):
    tree = Phylo.PhyloXMLIO.read(filename) #reads in phyloxml file as a Phylo tree
    tree0 = tree.phylogenies[0] 
    msa = tree0.to_alignment()
    
    f = open(os.path.splitext(filename)[0] + ".fasta", "w+")
    f.write(msa.format("fasta"))
    f.close()
