import numpy as np
import math
import re
import itertools
import pandas as pd
import time
import os
import pickle
import copy
from Functions import *
import cmcrameri.cm as cmc




def plot_data(dir_WAs, max_diff, min_N, start=0, stop=np.inf, step=1, **kwargs):
    N=start
    ls_names_met=['Method', 'Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps']
    
    while N<stop:
        Npath=os.path.join(dir_WAs,'N_'+str(N))
        if not os.path.exists(Npath):
            continue
        while diff<=max_diff:
            start_time = time.time()
            dpath=os.path.join(Npath,'diff_'+str(diff))
            if not os.path.exists(dpath):
                continue
            ls_met=[]
            full_methods=[]
            if diff==1:
                scrambler={1:np.arange(N)}
                WApath=os.path.join(dpath,'WAs')
                filenames = next(os.walk(WApath), (None, None, []))[2]
                for fname in filenames:
                    fdir=os.path.join(WApath,fname)
                    WA=np.genfromtxt(fdir, delimiter=",")
                    mean_exp, extra_exp,  _, perc_check= mean_metrics_precomp(well_assigner=WA,scrambler=scrambler, 
                                                                                differentiate=diff, **kwargs)
                    n_wells=WA.shape[1]
                    M_exp=np.round(mean_exp, 2)
                    max_comp=np.max(np.sum(WA, axis=0))
                    method=re.sub('^WA_', '', fname)
                    method=re.sub('_.*$', '', method)
                    ls_met.append([method, M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100])
                    full_methods.append(method)
                    
                Hier=calculate_metrics_hierarchical(n_compounds=N, differentiate=diff, **kwargs)
                ls_met.append(['Hierarchical']+ [np.round(i,2) for i in Hier[:-1]])
                full_methods.append('Hierarchical')
                df_met=pd.DataFrame(ls_met)
                #dft=pd.DataFrame({'Method':[a[0] for a in ls_met], 'Mean Experiments':[a[1] for a in ls_met],
                #                    'Max compunds':[a[2] for a in ls_met], 'N wells':[a[3] for a in ls_met],
                #                    'Percentage check':[a[4] for a in ls_met], 'Extra experiments':[a[5] for a in ls_met],
                #                    'Mean steps':[a[6] for a in ls_met],})



                idx_renamer={i:j for i,j in zip(df_met.index, full_methods)}
                col_renamer={i:j for i,j in zip(df_met.columns, ls_names_met)}
                df_met.rename(index=idx_renamer, columns=col_renamer, inplace=True)
                metname=os.path.join(dpath, 'Metrics_N_'+str(N)+'_diff_'+str(diff)+'.csv')
                #if not os.path.exists(dpath):
                #    os.makedirs(dpath)
                df_met.to_csv(metname)

                

                print('\n')
                print('-----------------------------------------------------')
                print("%s seconds required for N= %s and differentiate %s" % (np.round(time.time() - start_time, 1),N,diff))
                print('-----------------------------------------------------')
                diff+=1


        N+=step
                

