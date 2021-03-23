import pandas as pd
from os import listdir as ls
pd.set_option('display.max_columns', 50)
#brown sigsq line 7
#brown aicc line 24
#hansen alpha line 46
#hansen aicc line 68
#white noise sigsq line 79
#white noise aicc line 85 

def harvest_values(fstring):
    fn = open(direct+fstring,'r')
    lines = fn.readlines()
    fn.close()
    aiccs = []
    params = []
    sigsq_flag = True
    for line in lines:
        if '$aic.c' in line:
            aiccs.append(float(lines[lines.index(line)+1].split(' ')[-1]))
        elif 'AICc =' in line:
            aiccs.append(float(lines[lines.index(line)].split('\t')[-1].split(' ')[-1]))
        elif 'sigma.squared' in line:
            if sigsq_flag:
                params.append(float(lines[lines.index(line)+2].split(' ')[-1]))
                sigsq_flag=False
        elif '$alpha' in line:
            params.append(float(lines[lines.index(line)+2].split(' ')[-1]))
        elif 'sigsq' in line:
            params.append(float(lines[lines.index(line)].split('\t')[-1].split(' ')[-1]))
    if len(aiccs) == 3 and len(params)==3:
        return [aiccs, params]
    else:
        print([aiccs,params])
        return False

direct = 'AAC_uv_results_fixbl/'

residues = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','len']
colnames = ['linker']
for x in residues:
    colnames.append(x+'_best_fit')
    colnames.append(x+'_alpha')
    colnames.append(x+'_bm_sig')
    colnames.append(x+'_wn_sig')
log = open('test.log','w')
outdf = pd.DataFrame(columns=colnames)
for fstring in ls(direct):
        if fstring.split('_')[0] in residues:
                #fstring = 'A_model_fit_ENSGT00390000000002_1_Linker_1_2.txt'
                #log.write(fstring+'\n')
                fn = open(direct+fstring,'r')
                lines = fn.readlines()
                fn.close()

                if len(lines)>84: #only harvest results that finished for all 3 models, else an unfair comparison
                        residue = fstring.split('_')[0]
                        linker = fstring.split('.txt')[0].split('model_fit_')[1]
                        try:
                            vals = harvest_values(fstring)
                        except ValueError:
                            log.write(fstring)
                        if not vals:
                            log.write(fstring)
                        else:
                            aiccs, params = vals
                            alpha, bm_sig, wn_sig = params

                        '''
                        if residue !='len':
                            aiccs = [float(lines[23].split(' ')[-1]), float(lines[67].split(' ')[-1]), float(lines[84].split('\t')[-1].split(' ')[-1])]
                            alpha = float(lines[45].split(' ')[-1])
                            bm_sig = float(lines[6].split(' ')[-1])
                            wn_sig = float(lines[78].split('\t')[-1].split(' ')[-1])

                        else:
                            aiccs = [float(lines[23].split(' ')[-1]), float(lines[66].split(' ')[-1]), float(lines[94].split('\t')[-1].split(' ')[-1])]
                            alpha = float(lines[44].split(' ')[-1])
                            bm_sig = float(lines[6].split(' ')[-1])
                            wn_sig = float(lines[76].split('\t')[-1].split(' ')[-1])
                        '''

                        if min(aiccs)==aiccs[0]:
                            fit = 'brownian_motion'
                        elif min(aiccs)==aiccs[1]:
                            fit = 'ou'
                        elif min(aiccs)==aiccs[2]:
                            fit = 'white_noise'

                        
                        #entry = [linker, fit, float(lines[45].split(' ')[1]), float(lines[6].split(' ')[1]), float(lines[78].split('\t')[-1].split(' ')[-1])]
                        if linker not in outdf['linker'].tolist():
                            outdf = outdf.append({'linker': linker, residue+'_best_fit': fit, residue+'_alpha':alpha, residue+'_bm_sig':bm_sig, residue+'_wn_sig':wn_sig}, ignore_index=True)
                        else:
                            rownum = outdf.index[outdf['linker']==linker]
                            outdf.at[rownum,residue+'_best_fit'] = fit
                            outdf.at[rownum,residue+'_alpha'] = alpha
                            outdf.at[rownum,residue+'_bm_sig'] = bm_sig
                            outdf.at[rownum,residue+'_wn_sig'] = wn_sig


outdf.to_csv('AAC_model_params_all_residues_len_fixbl.csv', index=False)
