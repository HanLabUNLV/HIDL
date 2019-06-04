import os 
import simplejson as json
from sys import exit
import re

no_humans = []
f = open("domain_dict.json",'r')
DD = json.load(f)
direct = "/home/dbarth/Work/GeneTrees/SplitTrees/"
f.close()
for i in list(DD.keys()):
	try:
		parent = DD[i]["parent"][0].partition("_")[0]
	except IndexError:
		parent = DD[i]["parent"]	
	parent_list = []
	for j in [q for q in os.listdir(direct+parent) if q[-4:]==".fas"]:
		with open(direct+parent+"/"+j) as g:
			for data in g:
				if re.search(i,data):
					parent_list.append(j[:-10])
				elif re.search("ENSP0",data):
					print DD[i]['parent'] +" : "+ str(i)
					exit("human protein found, but not one we expected")
				else:
					no_humans.append(j[:-10])
	DD[i]["parent"]=parent_list
print "no human proteins found in the following gt"
print no_humans
f = open("domain_dict.json",'w')
f.write(json.dumps(DD))
f.close()
