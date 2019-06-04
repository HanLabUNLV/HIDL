import os
from sys import exit
from Bio import Phylo
from Bio import AlignIO
import json
import numpy as np
import re
from find_overlap import findOverlap

def align_to_numpy(alignment):
        array = np.array(list(alignment[0].seq))
        for rec in alignment[1:]:
                array = np.vstack([array, np.array(list(rec.seq))])
        return array

with open('shifted_dict.json','r') as f:
	dd = json.loads(f.readline())
def findGarbage():
	DD = {}
	garbage_proteins = []
	
	#starts off the dictionary with just the human data, this will be used as a guide to getting more data
	for i in list(dd.keys()): #i will be the species-specific protein_id
	        if len(dd[i]["start"])>1: #only grab the multi-domain proteins
	        	protein_dict = {}
	                d_len_matrix = np.subtract(np.array(dd[i]["end"]),np.array(dd[i]["start"]))
	                IAS_matrix = np.array([])
	                protein_dict["Domain Lengths"] = d_len_matrix
	                protein_dict["IAS Lengths"] = IAS_matrix
	
	                if dd[i]["parent"] not in DD: #some gt might have multiple human proteins
	                	DD[dd[i]["parent"]] = {}
	                        DD[dd[i]["parent"]][i] = protein_dict
	                else: #in the case of multiple human proteins, we don't want to reset the DD[i] dict
	                        DD[dd[i]["parent"]][i] = protein_dict	
	for i in list(DD.keys()): #i here is the gt id
		domain_lengths = []
		prot_keys = []
						
		for j in list(DD[i].keys()):	#first, get an array of the names of base sequences
			domain_lengths.append(DD[i][j]["Domain Lengths"])
			prot_keys.append(j)
		
		sorted_list = findOverlap(DD[i],dd)
		for j in range(len(sorted_list[0])):#have to add 0s to places where 1 protein has domain data and another in the same family doesn't
			for k in range(len(prot_keys)): #for each member of the base proteins
				if prot_keys[k] not in sorted_list[1][j]: #if member not in this domain (sorted_list[1] contains lists of proteins in each domain)
					domain_lengths[k] = np.insert(domain_lengths[k],j,0) #add a 0 to domain length, since it has none 	

		base_prot = [] #these will be the basis of further computations on alignments
		for j in range(len(domain_lengths[0])):
			min_dom = domain_lengths[0][j]
			min_prot = prot_keys[0]
			for k in range(len(domain_lengths)):
				if domain_lengths[k][j] < min_dom and domain_lengths[k][j] > 0:
					min_dom = domain_lengths[k][j]
					min_prot = prot_keys[k]	
			base_prot.append(min_prot)

		#second, bring in the whole alignment
	        alignment = AlignIO.read('SplitTrees/'+i.partition('_')[0]+'/' + i + ".fasta.fas",'fasta')
	        align_array = align_to_numpy(alignment)	

		#third, iterate over each of the alignments, starting from the base and calculating thereafter
		for k in range(len(dd[min_prot]['start'])): #for each domain (the kth domain)
			human_id = base_prot[k]
			start = dd[human_id]['start']
			end = dd[human_id]['end']
			for j in range(align_array.shape[0]): #j is the location of each protein id in the array
				print i + ' : ' + alignment[j].id
	                        protein_dict = {}
	                        protein_id = alignment[j].id #alignment[j].id should give the jth protein ID
	                        IAS_matrix = np.array([])
	                        d_len_matrix = np.array([])
	                        #in this block we calculate the IAS for the kth domain
	                        if k == 0: #in the first domain, there's no IAS on the left side
	                                IAS_number = 0
	                                Domain_number = DD[i][human_id]["Domain Lengths"][k]
	                                for l in range(end[k]-1,start[k+1]-1): #l runs over IAS characters
	                                        if align_array[j,l] != '-':
	                                                IAS_number += 1
	                                #walk back from end of the domain to give IAS penalty & Domain changes
	                                walk = end[k]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk!= start[k]-1): #deleted parts of domain are penalized
	                                        IAS_number -= 1
	                                        Domain_number -= 1
	                                        walk -= 1 #walk backwards from the end of the domain
	                                        ch = align_array[j,walk]
	                                #walk from the start of the domain to give ONLY domain changes
	                                walk = start[k]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk!= end[k]-1): #deleted parts of domain are penalized
	                                        Domain_number -= 1
	                                        walk += 1 #walk from the start of the domain
	                                        ch = align_array[j,walk]
	                                #save values
	                                IAS_matrix = np.append(IAS_matrix,[IAS_number])
	                                d_len_matrix = np.append(d_len_matrix,[Domain_number])
				
	                        if (k < len(start)-1) and (k>0): #only in the pairs between the first and last
	                                IAS_number = 0
					Domain_number = DD[i][human_id]["Domain Lengths"][k]
					#print (i,human_id)
	                                #calculate IAS score
	                                for l in range(end[k]-1,start[k+1]-1): #l runs over IAS characters
	                                        if align_array[j,l] != '-':
	                                                IAS_number += 1
	                                #walk from the start of the next domain to give IAS penalty on the right
	                                walk = start[k+1]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk<end[k]-1): #deleted parts of domain are penalized
	                                        IAS_number -= 1
	                                        walk += 1
	                                        ch = align_array[j,walk]
	
	                                #walk backwards from the end of this domain to penalize IAS *AND* change domain
	                                walk = end[k]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk !=start[k]-1):
        	                                IAS_number -= 1
	                                        Domain_number -= 1
	                                        walk -= 1
	                                        ch = align_array[j,walk]

                                	#walk from the start of this domain to change domain
                                	walk = start[k]-1
                                	ch = align_array[j,walk]
                                	while(ch == '-') and (walk !=end[k]):
                                	        Domain_number -= 1
                                	        walk += 1
                        	                ch = align_array[j,walk]
                	                #save values
        	                        IAS_matrix = np.append(IAS_matrix,[IAS_number])
	                                d_len_matrix = np.append(d_len_matrix,[Domain_number])
							
	                        if k == len(start)-1: #only in the last domain
	                                IAS_number = 0 #shouldn't ever get stored, used for debugging
	                                Domain_number = DD[i][human_id]["Domain Lengths"][k]
	                                #walk from the start of the domain to give Domain changes
	                                walk = start[k]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk!=end[k]-1): #deleted parts of domain are penalized
	                                        Domain_number -= 1
	                                        walk += 1
	                                        ch = align_array[j,walk]
	                                #walk back from the end of this domain to give changes
	                                walk = end[k]-1
	                                ch = align_array[j,walk]
	                                while (ch == '-') and (walk!=start[k]-1): #deleted parts of domain are penalized
	                                        Domain_number -= 1
	                                        walk -= 1
	                                        ch = align_array[j,walk]
	                                #save value
	                                d_len_matrix = np.append(d_len_matrix,[Domain_number])
				
	
				if Domain_number < ((0.5)*DD[i][human_id]["Domain Lengths"][k]):
					garbage_proteins.append((i,protein_id))
	return garbage_proteins

def delFastas(gt,protein_id):
	with open("/SplitTrees/"+gt.partition("_")[0]+'/' +gt+'.fasta','r+') as f: #needs to just be .fasta for final run
		string = '>'+protein_id
		data = f.readlines()
		f.seek(0)
		sequence_found = False
		for i in data:
			if sequence_found == True:
				if re.search('>',i)!=None:
					f.write(i)
					sequence_found=False
			elif re.search(string,i)!=None:
				sequence_found = True
			else:
				f.write(i)
		
		f.truncate()
		f.close()

def pruneNewick(gt,protein_id):
	tree = Phylo.read("/SplitTrees/"+gt.partition("_")[0]+'/' +gt+".newick",'newick')
	pruned_tree = tree.prune(target=protein_id)
	Phylo.write(tree,"/SplitTrees/"+gt.partition("_")[0]+'/' +gt+".newick",'newick')


need_cleaning = findGarbage()
for i in need_cleaning:
	delFastas(i[0],i[1])
	pruneNewick(i[0],i[1])







