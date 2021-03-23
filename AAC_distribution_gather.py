#This script will gather a distribution of each GAAC for each linker and pid
from os import listdir as ls
import pandas as pd
from multiprocessing import Pool

def read_pids(fn):
    pids = []
    for line in open(fn,'r'):
        pids.append(line.strip())
    return pids

linker_pid_dir = 'KCO_Linker_pids'
domain_pid_dir = 'KCO_Domain_pids'

linker_feature_dir = 'LinkerFastasFeatures'
domain_feature_dir = 'DomainFastasFeatures'

linker_length_dir = 'LinkerLengthsKCO'
domain_length_dir = 'DomainLengthsKCO'

masterdf = pd.DataFrame(columns=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y', 'pid','treeID','linker','domain','length'])

def process_file(ff):
    if 'Linker' in ff:
            full_linker = ff.split('.')[0]

            #initialize dataframe
            df = pd.read_csv(linker_feature_dir+'/'+ff, sep = '\t', header=None)
            df.drop(df.columns[0], axis=1, inplace=True) #drop rowname (0 for some mf reason on all of them)
            
            #create header for GAAC, add pids, ID, linker/domain
            df.columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            df['pid'] = read_pids(linker_pid_dir+'/'+full_linker+'_pids.txt')
            df['treeID'] = full_linker
            df['linker'] = True
            df['domain'] = False

            #initialize length column
            df['length'] = 0
            #get lengths for every pid and enter them 
            ldf = pd.read_csv(linker_length_dir + '/' +full_linker+'.lengths.txt', sep='\t')
            for index, row in df.iterrows():
                df.at[index,'length']=ldf.loc[ldf['seqId']==row['pid']]['totalNoGap']

            return df

    elif 'Domain' in ff:
            full_domain = ff.split('.')[0]

            #initialize dataframe
            df = pd.read_csv(domain_feature_dir+'/'+ff, sep = '\t', header=None)
            df.drop(df.columns[0], axis=1, inplace=True) #drop rowname (0 for some mf reason on all of them)
            
            #create header for GAAC, add pids, ID, linker/domain
            df.columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            df['pid'] = read_pids(domain_pid_dir+'/'+full_domain+'_pids.txt')
            df['treeID'] = full_domain
            df['linker'] = False
            df['domain'] = True

            #initialize length column
            df['length'] = 0
            #get lengths for every pid and enter them 
            ldf = pd.read_csv(domain_length_dir + '/' +full_domain+'.lengths.txt', sep='\t')
            for index, row in df.iterrows():
                df.at[index,'length']=ldf.loc[ldf['seqId']==row['pid']]['totalNoGap']

            return df

#add all linkers
linker_files = [x for x in ls(linker_feature_dir) if '.AAC.tsv' in x]
#add all domains
domain_files = [x for x in ls(domain_feature_dir) if '.AAC.tsv' in x]
#extend 
files = linker_files + domain_files

with Pool(75) as p:
    masterdflist = p.map(process_file, files)

masterdf = pd.concat(masterdflist, axis=0, ignore_index=True)
masterdf.to_csv('total_AAC_data.tsv',sep='\t',index=False)

