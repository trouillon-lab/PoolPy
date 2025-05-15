import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import pickle
import time
import os
import copy


from Functions import *


def is_consistent_precomp(well_assigner:np.array, differentiate:int, scrambler:dict) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    N=well_assigner.shape[0]
    for i in range(differentiate):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            this_sc=scrambler[diff]
            full_well_assigner=np.concatenate((full_well_assigner,np.bool(np.sum(well_assigner[this_sc], axis=1))))
    _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
    if len(counts)<full_well_assigner.shape[0]:
        return(False, full_well_assigner, counts)
    elif len(counts)==full_well_assigner.shape[0]:
        return(True,full_well_assigner, counts)
    else:
        print("Something is fishy")
        return(-1)
    
def is_consistent_precomp_alldiff(well_assigner:np.array, max_diff:int, scrambler:dict) -> list:
    if max_diff==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    N=well_assigner.shape[0]
    dict_counts={}
    for i in range(max_diff):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            this_sc=scrambler[diff]
            full_well_assigner=np.concatenate((full_well_assigner,np.bool(np.sum(well_assigner[this_sc], axis=1))))

        _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
        dict_counts.update({diff:counts})
    return dict_counts

def mean_metrics_precomp(well_assigner, differentiate, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent_precomp(well_assigner, differentiate, scrambler) 
    ET=extra_tests(counts)   
    rounds=np.sum(counts>1)/np.sum(counts>0)+1
    p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
    return BT+ET, ET,  rounds, p_check
    
#
def mean_metrics_precomp_alldiff(well_assigner, max_diff, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    dict_counts= is_consistent_precomp_alldiff(well_assigner, max_diff, scrambler) 
    dict_res={}
    for difo, counts in dict_counts.items():
        ET=extra_tests(counts)   
        rounds=np.sum(counts>1)/np.sum(counts>0)+1
        p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
        dict_res.update({difo:[BT+ET, ET,  rounds, p_check]})
    return dict_res


def sweep_metrics_precomp(dir_scramblers, dir_WAs, max_diff, start=50, stop=150, step=10, diff=2, **kwargs):
    N=start
    ls_names_met=['Method', 'Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps']
    while N<stop:
        Npath=os.path.join(dir_WAs,'N_'+str(N))
        diff=1
        if max_diff>1:
            del kwargs['differentiate']
            
            while diff<=max_diff:
                dpath=os.path.join(Npath,'diff_'+str(diff))
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
                    df_met.to_csv(metname)

                    diff+=1


                else:
                    print(diff)
                    this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                    scrambler.update({diff:this_sc_file})
                    dpath=os.path.join(Npath,'diff_'+str(diff),'WAs')
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
                    df_met.to_csv(metname)
                    
                    diff+=1

            N+=step
                    

        elif max_diff==1:
            N+=step
        else:
             N+=step



parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--dir_scramblers')
parser.add_argument('--dir_WAs')
parser.add_argument('--max_diff')
parser.add_argument('--max_dims')
parser.add_argument('--timeit')
parser.add_argument('--keep_ratios_constant')






args = parser.parse_args()



differentiate= 3 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
step= 10 if type(args.step)==type(None) else int(args.step)
keep_ratios_constant= False if type(args.keep_ratios_constant)==type(None) else args.keep_ratios_constant=='True' 
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'
max_diff= 4 if type(args.max_diff)==type(None) else int(args.max_diff)
max_dims= np.inf if type(args.max_dims)==type(None) else int(args.max_dims)



dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,'keep_ratios_constant': keep_ratios_constant,
             'start':start, 'stop':stop,  'step':step, 'dir_WAs':args.dir_WAs, 'dir_scramblers':args.dir_scramblers, 'max_diff': max_diff, 'max_dims':max_dims}


sweep_metrics_precomp(**dict_kwargs)
