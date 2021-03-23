from sklearn.decomposition import PCA
from sys import argv
from matplotlib import use
use('agg')
from matplotlib import pyplot as plt
import pandas as pd
from itertools import repeat as rep
import numpy as np 

def isLong(linker):
    if 'Linker' in linker:
        gtid = linker.split('_Linker')[0]
        linkerid = 'Linker' + linker.split('Linker')[-1]
        m = 
    else: 
        gtid = linker.split('_Domain')[0]
        linkerid = 'Domain' + linker.split('Domain')[-1]
    #print((gtid, linkerid, linker))


d_infile = argv[1]
l_infile = argv[2]
#outfile = 'feature_figures/'+'.'.join(d_infile.split('.')[1:-1]).split('/')[-1] + '.png'
outfile = 'feature_figures/'+'.'.join(d_infile.split('.')[1:-1]).split('/')[-1] + '_50_aa_min.png'
d_data = pd.read_csv(d_infile,sep='\t', header=None)
l_data = pd.read_csv(l_infile,sep='\t', header=None)
#print(d_data.shape)
dlabels = [x for x in rep('Domain',d_data.shape[0])]
llabels =  [x for x in rep('Linker',l_data.shape[0])]
labels = dlabels + llabels


data = pd.concat([d_data,l_data], axis=0)
print(data.head())

for index, row in data.iterrows():
    if isLong(row[0]):
        pass
    else:
        data.drop(index, inplace=True)

print(data.head())

exit()

pca=PCA(n_components=2)

#print(data.iloc[0:, 1:])
transformed = pca.fit_transform(data.iloc[0:, 1:-1])
plt.plot(transformed[d_data.shape[0]:,0], transformed[d_data.shape[0]:,1], label = 'Linker', alpha=0.2, marker = '.',linewidth=0, c='#f584f3')
plt.plot(transformed[0:d_data.shape[0],0], transformed[0:d_data.shape[0],1], label = 'Domain', alpha=0.2,marker='.',linewidth=0,c='#34c9eb')
plt.legend()
plt.title(' '.join(d_infile.split('.')[1:-1]).split('/')[-1] + ' PCA')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig(outfile)


