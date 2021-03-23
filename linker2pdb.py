

# read in protein IDs
dict_ensembl2pdb = {}
with open("ensemblID.linkers.pdb.ss.txt", "r") as file:
  for line in file:
    key, value = line.strip().split(" ")
    try:
      dict_ensembl2pdb[key].add(value)
    except KeyError:
      dict_ensembl2pdb[key] = set()
      dict_ensembl2pdb[key].add(value)
dict_tree = {}
dict_linker = {}
with open("protein2linkernss.txt", "r") as file:
  for line in file:
    prot, genetree, linker = line.strip().split(" ")
    try:
      dict_tree[genetree].add(linker)
    except KeyError:
      dict_tree[genetree] = set()
      dict_tree[genetree].add(linker)
    try:
      dict_linker[linker].add(prot)
    except KeyError:
      dict_linker[linker] = set()
      dict_linker[linker].add(prot)





for tree in dict_tree:
  print (tree)
  for linker in dict_tree[tree]: 
    print ("\t"+linker)
    for prot in dict_linker[linker]:
      print ("\t\t"+prot) 
      structs = dict_ensembl2pdb[prot]
      for pdb in structs:
        print ("\t\t\t"+pdb) 
        


