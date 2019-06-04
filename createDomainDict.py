import sys, os, json

famids = [x[0] for x in os.walk("geneDomains/")][1:] #the path to each gene tree folder
#print famids
domain_dict = {} #entries are formatted as {"human protein ID" : {'start':[start1,start2...], 'end':[end1,end2,...], 'parent': geneTreeID}}

for ii in famids:
	for filename in os.listdir(ii):
#		print ii +  filename
		with open(ii +"/"+ filename) as f:
			lst_o_dct = eval(f.readline()) #eval interprets datatypes, in this case a list of multiple dictionaries
			for dct in lst_o_dct:
				protein_id = dct["seq_region_name"] #search for the protein in the existing dictionary
				if domain_dict.get(protein_id) == None: #if a domain has not been recorded for this protein yet
					sub_dct = {}
					sub_dct["start"] = [dct["start"]]
					sub_dct["end"] = [dct["end"]]
					sub_dct["parent"] = ii[12:]
					domain_dict[protein_id] = sub_dct

				else: #if a domain has been recorded already for this protein
					domain_dict[protein_id]["start"].append(dct["start"])
					domain_dict[protein_id]["end"].append(dct["end"])
						
					domain_dict[protein_id]["start"].sort() #make sure the domains are in order
					domain_dict[protein_id]["end"].sort()
#	print domain_dict

with open("domain_dict.json", "w+") as d: #write the dictionary to a file
	d.write(json.dumps(domain_dict))
