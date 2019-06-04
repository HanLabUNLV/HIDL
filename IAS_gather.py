import numpy as np
import simplejson as json
import sys
from Bio import AlignIO
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend("agg")

#WHERE I LEFT OFF LAST:
#The first gene tree works fine
#The problem with the multiple human proteins was that the domains may not be the same between them
#I wrote some code to add gaps in the d_len_matrix (0s, they shouldn't get in the way [KW SHOULDNT])
#So wherever there's a domain that isn't covered by the other domains in a gt, this should fill those values in with 0
# therefore, I had to make sure that when choosing the smallest domain, that the domain lengths were > 0

#it works afaik, no further problems I can see
#i was wrong. i was very wrong.


def align_to_numpy(alignment):
        array = np.array(list(alignment[0].seq))
        for rec in alignment[1:]:
                array = np.vstack([array, np.array(list(rec.seq))])
        return array

class NumpyEncoder(json.JSONEncoder): #json dumps doesn't work with numpy, this should fix it
    def default(self, obj):		#there is no simpler way to put this dictionary into a human readable file format
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

dfile = open("shifted_dict.json", "r")
dd = eval(dfile.readline())
DD = {}

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
			
					
for i in ["ENSGT00530000063010"]: #list(DD.keys()):
	if len(list(DD[i].keys())) > 1:
		G = nx.Graph() #instantiate graph
		label_dict = {}
		for j in list(DD[i].keys()): #add nodes as tuples of start,finish
			for k in range(len(dd[j]["start"])):
				G.add_node((dd[j]["start"][k],dd[j]["end"][k]))
				label_dict[(dd[j]["start"][k],dd[j]["end"][k])] = j
		nx.set_node_attributes(G,"id",label_dict)
		Gnodes = list(G.nodes())
		for j in range(len(G.nodes())-1): #add edges if domains overlap
			for k in range(j,len(G.nodes())): #iterate through nodes not yet paired with j
				if Gnodes[j][0] > Gnodes[k][0] and Gnodes[j][0] < Gnodes[k][1]: #if the start of j is in range of k
					G.add_edge(Gnodes[j],Gnodes[k])
				elif Gnodes[j][1] > Gnodes[k][0] and Gnodes[j][1] < Gnodes[k][1]: #if the end of j is in range of k
					G.add_edge(Gnodes[j],Gnodes[k])
				elif Gnodes[k][0] > Gnodes[j][0] and Gnodes[k][0] < Gnodes[j][1]: #if the start of k is in range of j
					G.add_edge(Gnodes[j],Gnodes[k])	
				elif Gnodes[k][1] > Gnodes[j][0] and Gnodes[k][1] < Gnodes[j][1]: #if the end of k is in range of j
					G.add_edge(Gnodes[j],Gnodes[k])
		nx.draw(G,with_labels=True)
		plt.draw()
		plt.savefig("fig1.png")

		
		sorted_list = []
		for j in list(nx.connected_component_subgraphs(G)):
			if sorted_list == []:
				sorted_list.append(j.nodes())
			elif j.nodes()[0][0] < sorted_list[0][0][0]:
				sorted_list.insert(0,j.nodes())
			elif j.nodes()[0][0] > sorted_list[-1][0][0]:
				sorted_list.append(j.nodes())
			else:
				for k in range(len(sorted_list)):
					if j.nodes()[0][0] < sorted_list[k][0][0]:
						sorted_list.insert(k,j.nodes())
						break
		
		dom_num = 1
		for j in sorted_list: #j returns the list of tuples that cover a domain
			for k in list(DD[i].keys()):
				try:
					if (dd[k]["start"][dom_num],dd[k]["end"][dom_num]) not in j:
						print str((dd[k]["start"][dom_num],dd[k]["end"][dom_num])) +"not in "+str(j)
						DD[i][k]["Domain Lengths"]=np.insert(DD[i][k]["Domain Lengths"],dom_num,[0])			
				except IndexError:
					DD[i][k]["Domain Lengths"]=np.append(DD[i][k]["Domain Lengths"],[0])
			dom_num += 1
print sorted_list
for i in list(DD["ENSGT00530000063010"].keys()):
	print i + str(DD["ENSGT00530000063010"][i]["Domain Lengths"]) + str((dd[i]["start"],dd[i]["end"]))
sys.exit()	
DD_keys = list(DD.keys()) #needs to be instantiated here because more keys will be added in the for loop
for i in DD_keys: #i will be the gene-tree id
	print "Gene Tree ID : " + str(i)
	gt = i.partition("_")[0]
	human_id = list(DD[i].keys())[0] #returns the first key in the gene tree dictionary, there should be only human ids
	start = dd[human_id]["start"]
	end = dd[human_id]["end"]

	DD[i]["id_list"] = []
	DD[i]["IAS"] = np.array([])
	DD[i]["Domains"] = np.array([])
	
	aligndir = "/home/dbarth/Work/GeneTrees/SplitTrees/"
	alignment = AlignIO.read(aligndir + gt +'/' + i + ".fasta.fas",'fasta')
	align_array = align_to_numpy(alignment) 


	DD_i_keys = list(DD[i].keys())
	DD_i_keys.remove("Domains")
	DD_i_keys.remove("id_list")
	DD_i_keys.remove("IAS")

	for k in range(len(start)): #k runs over entries in start-end pairs
		#first, we pick the human protein with the smallest domain for this segment
		hid_score = DD[i][human_id]["Domain Lengths"][k]
                for entry in DD_i_keys:
	                entry_score = DD[i][entry]["Domain Lengths"][k]
			if entry_score < hid_score and entry_score > 0: 
                                human_id = entry
		#next we will iterate over all the proteins (remember this happens each domain)
		for j in range(align_array.shape[0]): #j is the location of each protein id in the array
			protein_dict = {} 
			protein_id = alignment[j].id[8:] #alignment[j].id[8:] should give the protein ID
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
				
				#calculate IAS score
				for l in range(end[k]-1,start[k+1]-1): #l runs over IAS characters
					if align_array[j,l] != '-': 
						IAS_number += 1
				#walk from the start of the next domain to give IAS penalty on the right
				walk = start[k+1]-1
				ch = align_array[j,walk]
				while (ch == '-') and (walk!=end[k]-1): #deleted parts of domain are penalized
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

			if protein_id not in DD[i]: #if the entry is not instantiated
				protein_dict["IAS Lengths"] = IAS_matrix
				protein_dict["Domain Lengths"] = d_len_matrix
				DD[i][protein_id] = protein_dict
			
			else: #if the entry *is* there
				if k < len(start)-1: #in all but the last domain, both IAS and domain are stored
					if protein_id[0:5] != "ENSP0": #if the protein is not a human one	
						new_ds = np.append(DD[i][protein_id]["Domain Lengths"],[Domain_number])
						DD[i][protein_id]["Domain Lengths"] = new_ds
	
						new_ias = np.append(DD[i][protein_id]["IAS Lengths"],[IAS_number])
						DD[i][protein_id]["IAS Lengths"] = new_ias
			
					else: #if it is a human protein, no need to re-apply domain information
						new_ias = np.append(DD[i][protein_id]["IAS Lengths"],[IAS_number])
						DD[i][protein_id]["IAS Lengths"] = new_ias
				else:	#in last domain, there's no IAS to the right, so no IAS stored
					if protein_id[0:5] != "ENSP0": #if the protein is not a human one	
						new_ds = np.append(DD[i][protein_id]["Domain Lengths"],[Domain_number])
						DD[i][protein_id]["Domain Lengths"] = new_ds	
			
	for protein_id in list(DD[i].keys()):
		if protein_id not in ("id_list","Domains","IAS"):
			DD[i]["id_list"].append(protein_id)
#			print("GT Extant Domain=" +str(DD[i]["Domains"].shape))
#			print("Adding Domains ... " + str(DD[i][protein_id]["IAS Lengths"]))
#			print(str(len(DD[i][protein_id]["Domain Lengths"])))
			#print protein_id + " " + dd[protein_id]["start"]
			if DD[i]["Domains"].shape != (0,):
				DD[i]["IAS"] = np.vstack((DD[i]["IAS"], DD[i][protein_id]["IAS Lengths"]))
				try:
					DD[i]["Domains"] = np.vstack((DD[i]["Domains"], DD[i][protein_id]["Domain Lengths"]))
				except ValueError:
					print DD[i]["Domains"]
					print DD[i][protein_id]["Domain Lengths"]
					print protein_id
					sys.exit()

			else:
				DD[i]["IAS"] = DD[i][protein_id]["IAS Lengths"]
				DD[i]["Domains"] = DD[i][protein_id]["Domain Lengths"]
	np.savetxt("/home/dbarth/Work/GeneTrees/IAS_Matrices/" + i + ".npy", DD[i]["IAS"])
	np.savetxt("/home/dbarth/Work/GeneTrees/DL_Matrices/" + i + ".npy", DD[i]["Domains"])

#	sys.exit("successful run through 1st gene tree")			
	
with open("master_dict.json", "w") as f:
	f.write(json.dumps(DD),cls = NumpyEncoder)



