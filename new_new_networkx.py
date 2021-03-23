from os import environ
environ['OPENBLAS_NUM_THREADS'] = '1'
import networkx as nx
from networkx.algorithms.approximation.clique import max_clique
from sys import exit, argv

#pass in gtid from command line (used for parallelizing)
gtid = argv[1]
#find data matching gtid
with open('master_dict_stratified.txt','r') as f:
    for line in f:
        if line.split('\t')[0]==gtid:
            dictline=line.strip()
            break
if dictline.split('\t')[1] == '{}':
    with open('clique_errors.log','a') as f:
        f.write('no entries for gt {}, exiting'.format(gtid)+'\n')
    exit(gtid)
index=0
#master = pickle.load(open('Pfam_sorted_shifted_dict_replacement.p','rb'))
master = {}
def readStratifiedDict(fn):
    #read in a stratified dict
    master = {}
    for line in fn:
        key, value = line.split('\t')
        master[key]=eval(value)
    return master
#master = readStratifiedDict(open('master_dict_stratified.txt','r'))
#test = 'ENSGT00390000012785_0'
out = open('master_dict_stratified_new.txt','a')
def cliqueToDomain(dictline): #for each gene tree
    dd, value = dictline.split('\t')
    master = {dd: eval(value)}
    G = nx.Graph()
    for ii in master[dd]: #for each sequence in this gene tree
        for jj in range(len(master[dd][ii])): #for each annotated domain datum
            if 'Linker' not in master[dd][ii][jj][0]:
                coordinates_and_ids = master[dd][ii][jj]
                G.add_node(ii+'.'+str(jj), coords=coordinates_and_ids[1:3], pid=coordinates_and_ids[3]) #add it to the graph, coords = start,end pair, pid=pfam id

    for kk in G.nodes(): 
        ens_id = kk.split('.')[0] #get the ID without the domain number

        start_kk = G.nodes[kk]['coords'][0] #start of the kth domain 
        end_kk = G.nodes[kk]['coords'][1] #end of the kth domain

        for mm in G.nodes(): 
            start_mm = G.nodes[mm]['coords'][0] #start of the mth domain
            end_mm = G.nodes[mm]['coords'][1] #end of the mth domain
            
            mid_kk = (start_kk+end_kk)/2
            mid_mm = (start_mm+end_mm)/2
            if (mid_kk<end_mm and mid_kk>start_mm) and (mid_mm<end_kk and mid_mm>start_kk):
                G.add_edge(mm,kk) #if their starts are closer than their ends 

    G_deep = G.copy()
    cliques = []
    #get list of nodes in cliques, remove those from the graph to solve problem of single domain member of two cliques
    while True:
        maxclique = list(max_clique(G))
        cliques.append(maxclique)
        G.remove_nodes_from(maxclique)
        if G.number_of_nodes() == 0:
            break

    domain_families = {} #Domain_0, Domain_1, etc. are keys, while values are a list of tuples like (ensembl_id, coords) 
    cliques.sort(key=lambda x: G_deep.nodes[x[0]]['coords'][0])
    #print([G_deep.nodes[x[0]]['coords'][0] for x in cliques])
    dom_idx = 0
    ddout = {}
    for i in cliques:
        for j in i:
            pid, domain = j.split('.')
            start_out, end_out = G_deep.nodes[j]['coords']
            pfam = G_deep.nodes[j]['pid']
            if pid in ddout:
                ddout[pid].append(['Domain_'+str(dom_idx), start_out, end_out, pfam]) 
            else:
                ddout[pid] = [['Domain_'+str(dom_idx), start_out, end_out, pfam]]
        dom_idx += 1
    return([dd,str(ddout)])
#print(dictline)
#print('\t'.join(cliqueToDomain(dictline))+'\n')
out.write('\t'.join(cliqueToDomain(dictline))+'\n')
'''
#parallelize the whole process, because the clique algorithm is NP hard        USE linux command parallel in py36 environment, cat the lines to a pipe into parallel 
pool=Pool(processes=85)
fin = open('master_dict_stratified.txt','r')
dictlines = fin.readlines()
for i in pool.map(cliqueToDomain, dictlines):
    print(i)    
    #out.write('\t'.join(i)+'\n')
'''

