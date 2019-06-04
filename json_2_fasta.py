# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:36:00 2018

@author: dylanbarth
"""

import json as j
import os

os.chdir("/DataDrives/") #changes directory to the DataDrives directory
#print os.getcwd() #prints current working directory as a check
direct = "/DataDrives/shared_107_dd3/clymer/structProject/familyIdJsons/" #where json files exist
new_files = "/home/dbarth/Work/familyIdFasta/" #place for new fasta files

for filename in os.listdir(direct): #iterates through the files in 'direct'
    os.chdir("/DataDrives/")
    with open(direct + filename) as jsonf:
        parsed = j.load(jsonf) #parse the json file using the built in python interpreter
        os.chdir(new_files) #move to new file location
        f = open(os.path.splitext(filename)[0] + ".fasta", "w+") #create new file
        uniprot_dict = parsed["MEMBERS"]["UNIPROT_proteins"] #parse for uniprot sequences
        for ii in uniprot_dict: #re-format into .fasta files
            f.write(">" + ii["protein_stable_id"] + "\n")
            f.write(ii["protein_alignment"] + "\n")
        ensembl_dict = parsed["MEMBERS"]["ENSEMBL_gene_members"] #parse for ensembl
        for jj in ensembl_dict: #re-format inot .fasta
            f.write(">" + ensembl_dict[jj][0]["protein_stable_id"] + "\n")
            f.write(ensembl_dict[jj][0]["protein_alignment"] + "\n")
        f.close()
