out = open('residue_model_preferences.txt','w')
aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','len']
models = ['BM','OU','WN']
props = {}
for aa in aas:
    bm = float(len(open(aa+'_data/'+aa+'_brownian_alphas.txt','r').readlines()))
    ou = float(len(open(aa+'_data/'+aa+'_ou_alphas.txt','r').readlines()))
    wn = float(len(open(aa+'_data/'+aa+'_wh_alphas.txt','r').readlines()))

    props[aa] = [str(bm/(bm+ou+wn)),str(ou/(bm+ou+wn)),str(wn/(bm+ou+wn))]

out.write('\t'.join(['Model','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','len\n']))
for x in (0,1,2):
    outlist = [models[x]]
    for aa in aas:
        outlist.append(props[aa][x])
    out.write('\t'.join(outlist)+'\n')

