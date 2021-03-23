fn = open('master_dict_stratified_newest.txt','r')
gt_num = 0
taxa = []
for line in fn:
    gt_num += 1
    gt, dd = line.split('\t')
    dd = eval(dd)
    linkers = {}
    for pid in dd:
        for entry in dd[pid]:
            if 'Linker_' in entry[0]:
                if entry[0] in linkers:
                    linkers[entry[0]] += 1
                else:
                    linkers[entry[0]] = 1
    for linker in linkers:
        taxa.append(linkers[linker])
out = open('num_taxa_in_linkers.txt','w')
out.write('taxa\n')
for i in taxa:
    out.write(str(i)+'\n')


