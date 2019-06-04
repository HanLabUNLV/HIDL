import numpy as np
from Bio import AlignIO
import json
from biopython_to_numpy import align_to_numpy
from sys import exit
direct = "/home/dbarth/Work/GeneTrees/SplitTrees/"

with open("domain_dict.json","r") as f: #un-shifted dictionary
	master_dict = eval(f.readline()) 
	progress = list(master_dict.keys())
	prog_count = 0.1
	for i in list(master_dict.keys()): #i is the human protein associated with each gene tree
		gt = master_dict[i]["parent"][0] 
		gt = gt.partition("_")[0] #gt is the gene tree id without the subscript we placed on there
		print gt +" : " + i
		print direct + gt +"/" + master_dict[i]["parent"][0] + ".fasta.fas"
		alignment = AlignIO.read(direct + gt +"/" + master_dict[i]["parent"][0] + ".fasta.fas",'fasta')
		align_array = align_to_numpy(alignment)
			
		a = 0
		for record in alignment:
			if record.id == i:
				human_index = a
				break	
			a += 1
		
		walk = 0 
		dashes = 0
		s_index = 0
		e_index = 0
		start_add = []
		end_add = []
							
		while walk <= (master_dict[i]["end"][-1]-1):
			if align_array[human_index,walk] == '-':
				dashes += 1
			if (s_index==e_index) and  (walk >= (master_dict[i]["start"][s_index]-1)):
				start_add.append(dashes)
				s_index += 1
			if (walk >= (master_dict[i]["end"][e_index]-1)) and (s_index >= e_index):
				end_add.append(dashes)
				e_index += 1	
			walk += 1


		for j in range(len(start_add)):
			master_dict[i]["start"][j] += start_add[j]
			master_dict[i]["end"][j] += end_add[j]

		progress.remove(i)
		if round(1 - len(progress)/len(list(master_dict.keys())), 1) > prog_count:
			print "progress: " + str(prog_count*100) + "%"
			prog_count += 0.1

with open("shifted_dict.json","w") as g:
	g.write(json.dumps(master_dict))

