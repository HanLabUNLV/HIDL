from Bio import Phylo
import simplejson as json
import re
import os
import gc
from sys import exit
import re


direct = "/home/dbarth/Work/GeneTrees/geneJsons/" #directory where full gt jsons exist
splt_dir = "/home/dbarth/Work/GeneTrees/SplitTrees/" #directory where new data will be stored
os.chdir(direct)
gt_ids = [] #gene trees with nodes who are older than mya_limit 
mya_limit = 400

def bfs(sub_tree): #simple breadth first search of the .json gene tree
	if ('timetree_mya' in sub_tree["taxonomy"].keys()) and (sub_tree["taxonomy"]["timetree_mya"] > mya_limit): #if this node is still too old
		if 'children' in list(sub_tree.keys()):
			for sub in sub_tree["children"]: #try all this nodes children nodes
				bfs(sub)		 
		else: #except if this node has no children
			pass
	elif 'timetree_mya' not in sub_tree["taxonomy"].keys():	
		if 'children' in list(sub_tree.keys()):
			for sub in sub_tree["children"]: #try all this nodes children nodes
				bfs(sub)		 
		else: #except if this node has no children
			pass
	else: #exit condition met, sub_tree is younger than mya_limit
		trees.append(sub_tree)		



for filename in os.listdir(direct): #iterate through the json files in direct
	f = open(filename, 'r') 
	data = f.readline()
	if data[0] != "{": #sometimes ensembl puts the gt id as the first line, so just read the next line if they did
		data = f.readline()
	r = re.findall(r'timetree_mya\': (.*?)\}',data) #start with the timetree keys and its values
	r = [i for i in r if i.strip('0123456789') == '.'] #grab the number from the values
	r = [float(i) for i in r] #turn that number into a float
	r = [i for i in r if i > mya_limit] #only keep trees that are older than mya_limit, otherwise what's the point of splitting them?
	
	if r != []:
		try:
			os.mkdir(splt_dir + filename[:-5] +"/") #make the directory for this gene tree if it doesn't exist
		except OSError:
			pass
		gt_ids.append(filename[:-5]) #gt w nodes - '.json' file extension
	elif re.search('Human',data) != None : #if this gt's oldest node was younger than mya_limit, and has human data 
		try:
			os.mkdir(splt_dir + filename[:-5] +"/") #make the directory for this gene tree if it doesn't exist
		except OSError:
			pass

                data = re.sub('u\"(.+?)\'(.+?)\",','\"\",',data)  #use regex to clean text bc json is picky
                data = re.sub('u\'','\"',data)                    #^
                data = re.sub('u\"','\"', data)                   #^
                data = re.sub('\'','\"',data)   
		
		data = json.loads(data)

		g = open(splt_dir + filename[:-5] +"/"+filename,'w') #just write the file
		g.write(json.dumps(data))
		g.close()


	f.close()
for filename in gt_ids: #for our filtered gene trees
	trees = [] #this will hold all the sub-trees
	with open(filename+'.json','r') as f:
		print filename #for debugging and progress 
		dic = f.readline() 
		if dic[0] != "{": #if the first line wasn't a dictionary (ensembl format)
			dic = f.readline()
		else:
			f.seek(0)
			dic = f.readline()	
		dic = re.sub('u\"(.+?)\'(.+?)\",','\"\",',dic)	#use regex to clean text bc json is picky
		dic = re.sub('u\'','\"',dic) 			#^
		dic = re.sub('u\"','\"', dic) 			#^
		dic = re.sub('\'','\"',dic) 			#^
#		print dic[47058:47258]

		dic = json.loads(dic) #load our sanitized string as a dictionary	

		bfs(dic["tree"]) #call the breadth-first search, which will append sub-trees to trees		
	
		h_trees = [] #trees with human seq in them
		if trees != []: 
			for i in trees: #iterate through all sub-trees	
				if re.search('Human',str(i)) != None: #if the sub-tree doesn't have human protein, we discard
					h_trees.append(i)
		
		for i in range(len(h_trees)): #write the dic to files
			g = open(splt_dir+filename+"/"+filename + "_" + str(i)+".json",'w')
			g.write(json.dumps(h_trees[i]))
			g.close()
	 

#things to do:
#	1. prune sub trees breadth first
# 	2. save these trees as a json file
#	3. use the timetree data for creating newicks
#	4. extract fasta from these new jsons
