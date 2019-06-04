import simplejson as json
import os
from sys import exit
def json_2_newick(json_dict):
	if 'children' not in list(json_dict.keys()):	
		try:
			if isinstance(json_dict["id"],list):
				spec_id = json_dict["id"][0]["accession"]
			elif isinstance(json_dict["id"],dict):	
				spec_id = json_dict["id"]["accession"]
			elif isinstance(json_dict["id"],str):
				spec_id = json_dict["id"]
			else:
				print "id key returns unrecognized data type"
				print json_dict["id"]
				exit()
		except KeyError:
			spec_id = json_dict["taxonomy"]["scientific_name"]	
		
		try:
			branch_len = json_dict["branch_length"]
		except KeyError:
			branch_len = 0.		
#		print spec_id + ':' + str(branch_len)
		return spec_id +':'+ str(branch_len)
	else:
		
		try:
			branch_len = json_dict["branch_length"]
		except KeyError:
			branch_len = 0.	
		try:
			if isinstance(json_dict["id"],list):
				spec_id = json_dict["id"][0]["accession"]
			elif isinstance(json_dict["id"],dict):	
				spec_id = json_dict["id"]["accession"]
			elif isinstance(json_dict["id"],str):
				spec_id = json_dict["id"]
			else:
				print "id key returns unrecognized data type"
				print json_dict["id"]
				exit()
		except KeyError:
			spec_id = json_dict["taxonomy"]["scientific_name"]	

#		print spec_id + ':' + str(branch_len)
		children = []
		for sub in json_dict['children']:
			new = json_2_newick(sub)
			if isinstance(new,list):
				new = '['+','.join(new)+']'
			children.append(new)
		chstr = '['+','.join(children)+']'
		return [spec_id, '['+','.join(children)+']'+':'+str(branch_len)]
				

direct = "/home/dbarth/Work/GeneTrees/SplitTrees/"
for d in os.listdir(direct):
	for filename in os.listdir(direct + d):
		if filename[-5:] == ".json":
			with open(direct + d +"/"+ filename) as g:
				data = json.load(g)
				try:
					data =  data["tree"]
				except KeyError:
					pass
				newick = json_2_newick(data)
				newick = str(newick) + ';'
				newick = newick.replace('[','(')
				newick = newick.replace(']',')')
				newick = newick.replace('\'','')
			with open(direct + d +"/"+ filename[:-5] + ".newick",'w') as f:
				f.write(str(newick))


#json_2_newick seems to work
#TODO: iterate through a directory, find all directories, then iterate through files in those directories
