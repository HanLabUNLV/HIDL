import numpy as np
import simplejson as json
import sys
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend("agg")

#this program takes 2 dictionaries DD and dd, the former with structure {gt_id:prot_id{}....}, ONLY 1 ENTRY / Gene Tree
#the latter with structure {prot_id: {start:[],end:[],parent:gt_id}....}
#i'm not yet sure this returns the right things... or inputs the right things....
def findOverlap(DD_i,dd):
	if len(list(DD_i.keys())) > 1:
		G = nx.Graph() #instantiate graph
		label_dict = {}
		for j in list(DD_i.keys()): #add nodes as tuples of start,finish
			for k in range(len(dd[j]["start"])):
				G.add_node((dd[j]["start"][k],dd[j]["end"][k]))
				if (dd[j]["start"][k],dd[j]["end"][k]) not in label_dict.keys():
					label_dict[(dd[j]["start"][k],dd[j]["end"][k])] = [j]
				else:
					label_dict[(dd[j]["start"][k],dd[j]["end"][k])].append(j)
		nx.set_node_attributes(G,"id",label_dict)
		Gnodes = list(G.nodes())
		for j in range(len(G.nodes())-1): #add edges if domains overlap
			for k in range(j,len(G.nodes())): #iterate through nodes not yet paired with j
				if Gnodes[j][0] >= Gnodes[k][0] and Gnodes[j][0] <= Gnodes[k][1]: #if the start of j is in range of k
					if abs(Gnodes[j][0] - Gnodes[k][0]) < abs(Gnodes[j][0] - Gnodes[k][1]): #if jstart is closer to kstart than kend
						G.add_edge(Gnodes[j],Gnodes[k])
				elif Gnodes[j][1] >= Gnodes[k][0] and Gnodes[j][1] <= Gnodes[k][1]: #if the end of j is in range of k
					if abs(Gnodes[j][1] - Gnodes[k][1]) < abs(Gnodes[j][1] - Gnodes[k][0]): #if jend is closer to kend than kstart
						G.add_edge(Gnodes[j],Gnodes[k])
				elif Gnodes[k][0] >= Gnodes[j][0] and Gnodes[k][0] <= Gnodes[j][1]: #if the start of k is in range of j
					if abs(Gnodes[j][0] - Gnodes[k][0]) < abs(Gnodes[j][0] - Gnodes[k][1]): #if kstart is closer to jstart than jend
						G.add_edge(Gnodes[j],Gnodes[k])	
				elif Gnodes[k][1] >= Gnodes[j][0] and Gnodes[k][1] <= Gnodes[j][1]: #if the end of k is in range of j
					if abs(Gnodes[j][1] - Gnodes[k][1]) < abs(Gnodes[j][1] - Gnodes[k][0]): #if kend is closer to jend than jstart
						G.add_edge(Gnodes[j],Gnodes[k])
#		nx.draw(G,with_labels=True)
#		plt.draw()
#		plt.savefig("fig1.png")
	
			
		sorted_list = []
		labels = []
		print "# of subgraphs" , len(list(nx.connected_component_subgraphs(G)))
		for j in list(nx.connected_component_subgraphs(G)):
			if sorted_list == []:
				sorted_list.append(j.nodes())
				labels.append([label for sublist in nx.get_node_attributes(j,'id').values() for label in sublist])
			elif j.nodes()[0][0] < sorted_list[0][0][0]:
				sorted_list.insert(0,j.nodes())
                                labels.insert(0,[label for sublist in nx.get_node_attributes(j,'id').values() for label in sublist])
			elif j.nodes()[0][0] > sorted_list[-1][0][0]:
				sorted_list.append(j.nodes())
                                labels.append([label for sublist in nx.get_node_attributes(j,'id').values() for label in sublist])
			else:
				for k in range(len(sorted_list)):
					if j.nodes()[0][0] < sorted_list[k][0][0]:
						sorted_list.insert(k,j.nodes())
						labels.insert(k,[label for sublist in nx.get_node_attributes(j,'id').values() for label in sublist])
						break
			
#		dom_num = 1
#		for j in sorted_list: #j returns the list of tuples that cover a domain
#			for k in list(DD[i].keys()):
#				try:
#					if (dd[k]["start"][dom_num],dd[k]["end"][dom_num]) not in j:
#						print str((dd[k]["start"][dom_num],dd[k]["end"][dom_num])) +"not in "+str(j)
#						DD[i][k]["Domain Lengths"]=np.insert(DD[i][k]["Domain Lengths"],dom_num,[0])			
#				except IndexError:
#					DD[i][k]["Domain Lengths"]=np.append(DD[i][k]["Domain Lengths"],[0])
#			dom_num += 1
	else:
		sorted_list = []
		labels = []

	return [sorted_list,labels]
