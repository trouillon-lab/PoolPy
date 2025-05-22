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



def get_max_C(n_compounds, max_compounds):
    return int(n_compounds/2) if max_compounds==0 else max_compounds
def get_min_C(n_compounds, MC):
    return int(np.sqrt(n_compounds)) if int(np.sqrt(n_compounds))<MC else int(MC/2)

def get_max_W(n_compounds):
    return int(np.log2(n_compounds))
def get_min_W(n_compounds):
    return int(2*np.sqrt(n_compounds))


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


def mean_metrics_precomp(well_assigner, differentiate, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent_precomp(well_assigner, differentiate, scrambler) 
    ET=extra_tests(counts)   
    rounds=np.sum(counts>1)/np.sum(counts>0)+1
    p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
    return BT+ET, ET,  rounds, p_check



def find_rand_params_precomp(n_compounds:int, n_compounds_per_well=0, n_wells=0, guesses=0, 
                     max_compounds=0, max_redundancy=4, min_redundancy=1,**kwargs):
    skip_compounds=True
    skip_wells=True
    if n_compounds_per_well==0:
        skip_compounds=False
    if n_wells==0:
        skip_wells=False
    if guesses==0:
        guesses=n_compounds

    MC= get_max_C(n_compounds, max_compounds)
    mc=get_min_C(n_compounds, MC)
    arr_comp=np.arange(int(mc),int(MC+1))
    mw=get_min_W(n_compounds)
    MW=get_max_W(n_compounds)
    while MW-mw<10:
        mw=int(abs(mw-1))
        MW=int(MW+1)

    arr_wells=np.arange(mw,MW)
    min_tests=np.inf
    N_tries=0
    for comp in arr_comp:
        if skip_compounds:
            comp=n_compounds_per_well

        for wells in arr_wells:
            if skip_wells:
                if skip_compounds:
                    
                    return 0
                wells=n_wells
                
            if comp*wells>max_redundancy*n_compounds*np.log2(n_compounds) or comp*wells<min_redundancy*n_compounds: continue 
            WA_tmp, mean_exp, p_check=evaluate_rand_design(n_compounds=n_compounds,n_wells=wells,n_compounds_per_well=comp, return_me=True, guesses=guesses, **kwargs)
            N_tries+=1
            if mean_exp<min_tests:
                Comp=comp
                Wells=wells
                min_tests=mean_exp
                min_wa=WA_tmp
                min_pcheck=p_check
            if skip_wells:
                break
        if skip_compounds:
            break

    print('\n')
    print('----------------------------------------------------------------------------------------------------------')
    print("Evaluated %s different random designs each with %s configurations " % (N_tries, guesses+1))
    print('\n')

    return Comp, Wells, min_tests, min_wa, min_pcheck


def evaluate_rand_design(n_compounds:int,  differentiate:int,scrambler:dict, n_compounds_per_well=0, 
                        n_wells=0, guesses=0,  return_me=False, **kwargs):
    min_tests=np.inf
    second_axis=np.tile(np.arange(n_wells),n_compounds_per_well).reshape(n_compounds_per_well,-1)
    for i in range(guesses):
        #print(i)
        idt=np.random.randint(0,n_compounds,size=(n_compounds_per_well,n_wells) )
        well_assigner=np.zeros((n_compounds,n_wells))==1
        well_assigner[idt, second_axis]=True
        if guesses==1:
            if return_me:
                mean_exp, _, _, p_check= mean_metrics_precomp(well_assigner=well_assigner, 
                                                            differentiate=differentiate,scrambler=scrambler,**kwargs)
                return well_assigner, mean_exp, p_check
            return well_assigner
        mean_exp, _, _, p_check= mean_metrics_precomp(well_assigner=well_assigner,
                                                    differentiate=differentiate, scrambler=scrambler, **kwargs)
        if p_check<1:
            if return_me:
                return well_assigner,  mean_exp, p_check
            return well_assigner
        elif mean_exp<min_tests: 
            best_wa=well_assigner.copy()
            min_tests=mean_exp
            min_pcheck=p_check

    if return_me:
        return best_wa,  min_tests, min_pcheck
    
    return best_wa


def assign_wells_random_precomp(n_compounds:int,  differentiate:int,scrambler:dict, n_compounds_per_well=0, 
                        n_wells=0, guesses=0, Evaluate=False, return_me=False, **kwargs)->np.array:

    _,_, min_tests, WA_rand, p_check=find_rand_params_precomp(n_compounds=n_compounds, differentiate=differentiate, 
                                 n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, 
                                 guesses=guesses,scrambler=scrambler, **kwargs)
    if return_me:
        return WA_rand,  min_tests, p_check
    
    return WA_rand




def rand_sweep_diff(n_compounds, max_diff, dir_scramblers, Npath, **kwargs):
    if 'differentiate' in kwargs.keys():
        del kwargs['differentiate']
    N=n_compounds
    
    if max_diff>1:

        for di in range(max_diff):
            start_time = time.time()
            diff=di+1
            if diff==1:
                dpath=os.path.join(Npath,'diff_'+str(diff))
                WApath=os.path.join(dpath,'WAs')
                scrambler={1:np.arange(n_compounds)}
                WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                                differentiate=diff,scrambler=scrambler, return_me=True, **kwargs )
                extra_exp=WA_rand.shape[1]+min_tests
                #.append(['Random', min_tests, np.max(np.sum(WA_rand, axis=0)), WA_rand.shape[0], int(perc_check),  extra_exp,1+perc_check/100])
                full_file_dir=os.path.join(dpath,'Random_diff_'+str(diff)+'_NS_'+
                                               str(n_compounds)+'_NW_'+str(WA_rand.shape[0])+
                                               '_MS_'+str(np.max(np.sum(WA_rand, axis=0)))+
                                                '_PC_'+ str(int(perc_check)) +'_EE_'+str(extra_exp)+".txt")
                if not os.path.exists(dpath):
                    os.makedirs(dpath)
                open(full_file_dir, 'a').close()
                if not os.path.exists(WApath):
                    os.makedirs(WApath)
                thisfile=os.path.join(WApath,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+'.csv')
                np.savetxt(thisfile, WA_rand.astype(bool), delimiter=",")


            else:
                dpath=os.path.join(Npath,'diff_'+str(diff))
                WApath=os.path.join(dpath,'WAs')
                this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                this_scrambler=np.load(this_sc_file)['sc']
                scrambler.update({diff:this_scrambler})
                WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                                differentiate=diff,scrambler=scrambler, return_me=True, **kwargs  )
                extra_exp=WA_rand.shape[1]+min_tests
                #.append(['Random', min_tests, np.max(np.sum(WA_rand, axis=0)), WA_rand.shape[0], int(perc_check),  extra_exp,1+perc_check/100])
                full_file_dir=os.path.join(dpath,'Random_diff_'+str(diff)+'_NS_'+
                                               str(n_compounds)+'_NW_'+str(WA_rand.shape[0])+
                                               '_MS_'+str(np.max(np.sum(WA_rand, axis=0)))+
                                                '_PC_'+ str(int(perc_check)) +'_EE_'+str(extra_exp)+".txt")
                if not os.path.exists(dpath):
                    os.makedirs(dpath)
                open(full_file_dir, 'a').close()
                if not os.path.exists(WApath):
                    os.makedirs(WApath)
                thisfile=os.path.join(WApath,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+'.csv')
                np.savetxt(thisfile, WA_rand.astype(bool), delimiter=",")

            print("%s seconds required for N= %s and differentiate %s" % (np.round(time.time() - start_time, 1),n_compounds,diff))
            print('----------------------------------------------------------------------------------------------------------')



def rand_N_sweep(start, stop, step,dir_WAs, **kwargs):
    n_compounds=start
    while n_compounds<=stop:
        Npath=os.path.join(dir_WAs,'N_'+str(n_compounds))
        rand_sweep_diff(n_compounds=n_compounds, Npath=Npath, **kwargs)
        n_compounds+=step






parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--dir_scramblers')
parser.add_argument('--dir_WAs')
parser.add_argument('--max_diff')
parser.add_argument('--timeit')
parser.add_argument('--max_compounds')
parser.add_argument('--n_compounds_per_well')
parser.add_argument('--n_wells')
parser.add_argument('--rand_guesses')







args = parser.parse_args()




differentiate= 3 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
step= 10 if type(args.step)==type(None) else int(args.step)
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'
max_diff= 4 if type(args.max_diff)==type(None) else int(args.max_diff)
rand_guesses= 10 if type(args.rand_guesses)==type(None) else int(args.rand_guesses)



dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,'guesses':rand_guesses,
             'start':start, 'stop':stop,  'step':step, 'dir_WAs':args.dir_WAs, 'dir_scramblers':args.dir_scramblers, 'max_diff': max_diff,}

if type(args.max_compounds)!=type(None): 
    dict_kwargs.update({'max_compounds':int(args.max_compounds)})
if type(args.n_compounds_per_well)!=type(None): 
    dict_kwargs.update({'n_compounds_per_well':int(args.n_compounds_per_well)})
if type(args.n_wells)!=type(None): 
    dict_kwargs.update({'n_wells':int(args.n_wells)})

rand_N_sweep(**dict_kwargs)


