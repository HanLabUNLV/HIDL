from Bio import Phylo

def get_parent(tree,child): #returns one node above the node given
        node_path = tree.get_path(child)
        return node_path[-2]

def get_children(tree,node): #returns list of all direct descendent nodes
        children = []
        node_clade = tree.find_clades(target=node)
        for i in node_clade:
                for j in i:
                        children.append(j.name)
        return children

def is_sibling(tree,child1,child2): #returns True if child1 and child2 are siblings
        btwn = tree.trace(child1,child2) #trace is supposed to tell me all the nodes between the values excluding the start and finish
        for i in btwn:                  #however, sometimes it includes start and finish, so... we gotta remove em
                if (i.name == child1) or (i.name == child2):
                        btwn.remove(i)
        if len(btwn)==1:
                return True
        else:
                return False

def generate_ids(tree): #some internal nodes do not have a name attribute, but do have a unique identifier, we fix it with this
        for idx, clade in enumerate(tree.find_clades()): #NOTE: internal nodes that are renamed are not always unique to the tree
                if not clade.name:
                        clade.name = "node_"+str(idx)

def get_siblings(tree,node): #returns list of siblings -- direct descendants of parent
        sibs = []
        parent = get_parent(tree,node).name
        for i in get_children(tree,parent):
                if i!=node:
                        sibs.append(i)
        return sibs


