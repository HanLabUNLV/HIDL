import pandas as pd
from statistics import mean
infile = 'AAC_model_params.csv'
linker_length_dir = 'LinkerLengthsKCO'

df = pd.read_csv(infile)

for index, row in df.iterrows():
    linker = row.loc['linker']
    ldf = pd.read_csv(linker_length_dir + '/' +linker+'.lengths.txt', sep='\t')
    df.at[index,'length']=mean(ldf['totalNoGap'])

df.to_csv('AAC_model_params_lens.csv', index=False)
