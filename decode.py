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
import json


from Functions import *


def decode_precomp(well_assigner:np.array, differentiate:int, 
                   scrambler:dict, readout:np.ndarray, max_differentiate=-1, sweep=False, **kwargs) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    if max_differentiate<1:
        N=well_assigner.shape[0]
        sc_list=np.arange(N).tolist()
        for i in range(differentiate):
            diff=i+1
            if diff ==1:
                full_well_assigner=well_assigner.copy()
            else:
                this_sc=scrambler[diff]
                #print(this_sc)
                #print(well_assigner)
                #print(diff)
                full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
                sc_list.extend(this_sc.tolist)
        #outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=True)
        
        if sweep:
            outcome_dict={}
            outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=False).astype(bool)
            for outcome in outcomes:
                idxs = np.all(outcome == full_well_assigner, axis=1)
                outcome_dict.update({tuple_to_str(tuple(outcome)):list(itertools.compress(sc_list,idxs))})
                print(sc_list)
                print(idxs)
                print(outcome)
                print(outcome_dict)
                
            exit()
                
            return outcome_dict

        else:
            idxs = np.all(outcome == full_well_assigner, axis=1)
            return list(itertools.compress(sc_list,idxs))
        
    else:
        full_od={}
        N=well_assigner.shape[0]
        sc_list=[np.arange(N)]
        for differentiate in range(max_differentiate):

            diff=differentiate+1
            if diff ==1:
                full_well_assigner=well_assigner.copy()
            else:
                this_sc=scrambler[diff]
                #print(this_sc)
                #print(well_assigner)
                #print(diff)
                full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
                sc_list.extend(this_sc.tolist)
        #outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=True)
        
            if sweep:
                outcome_dict={}
                outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=False).astype(bool)
                for outcome in outcomes:
                    idxs=np.prod(outcome==full_well_assigner, axis=1)
                    outcome_dict.update({tuple_to_str(tuple(outcome)):list(itertools.compress(sc_list,idxs))})
                full_od.update({diff:outcome_dict})

            else:
                idxs=np.prod(readout==full_well_assigner, axis=1)
            full_od.update({diff:list(itertools.compress(sc_list,idxs))})
        return full_od
    
def decode_sweep(dir_scramblers, dir_WAs, differentiate:int,
            max_differentiate,
            start=50, stop=150, step=10,
             **kwargs) -> list:
    N=start
    while N<stop:
        Npath=os.path.join(dir_WAs,'N_'+str(N))
        diff=1
        if max_differentiate>=1:

            diff=1
            if 'differentiate' in kwargs.keys():
                del kwargs['differentiate']
            while diff<=max_differentiate:
                start_time = time.time()
                dpath=os.path.join(Npath,'diff_'+str(diff))
                if diff==1:
                    scrambler={1:np.arange(N)}
                    WApath=os.path.join(dpath,'WAs')
                    filenames = next(os.walk(WApath), (None, None, []))[2]
                    for fname in filenames:
                        fdir=os.path.join(WApath,fname)
                        WA=np.genfromtxt(fdir, delimiter=",")
                        method=re.sub('^WA_', '', fname)
                        method=re.sub('_.*$', '', method)
                        dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 
                        scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                        decpath=os.path.join(dpath,'decoders')
                        decname=os.path.join(decpath, 'decoder_'+method+'.json')
                        print(dict_decode)
                        json.dump( dict_decode, open(decname, 'w' ) )

                else:
                    this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                    this_scrambler=np.load(this_sc_file)['sc']
                    scrambler.update({diff:this_scrambler})
                    dpath=os.path.join(Npath,'diff_'+str(diff))
                    WApath=os.path.join(dpath,'WAs')
                    filenames = next(os.walk(WApath), (None, None, []))[2]
                    for fname in filenames:
                        fdir=os.path.join(WApath,fname)
                        WA=np.genfromtxt(fdir, delimiter=",")
                        method=re.sub('^WA_', '', fname)
                        method=re.sub('_.*$', '', method)
                        dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 

                        scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                        decpath=os.path.join(dpath,'decoders')
                        decname=os.path.join(decpath, 'decoder_'+method+'.json')

                        json.dump( dict_decode, open(decname, 'w' ) )
                diff+=1

        else:
            diff=1    
            while diff<=differentiate:  
                if diff==1:
                        scrambler={1:np.arange(N)}
                        WApath=os.path.join(dpath,'WAs')
                        filenames = next(os.walk(WApath), (None, None, []))[2]
                        for fname in filenames:
                            fdir=os.path.join(WApath,fname)
                            WA=np.genfromtxt(fdir, delimiter=",")
                            method=re.sub('^WA_', '', fname)
                            method=re.sub('_.*$', '', method)
                            dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 
                            scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                            decname=os.path.join(dpath, 'decoder_'+method+'.json')
                            json.dump( dict_decode, open(decname, 'w' ) )

                else:
                    this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                    this_scrambler=np.load(this_sc_file)['sc']
                    scrambler.update({diff:this_scrambler})
                    dpath=os.path.join(Npath,'diff_'+str(diff))
                    WApath=os.path.join(dpath,'WAs')
                    filenames = next(os.walk(WApath), (None, None, []))[2]
                    for fname in filenames:
                        fdir=os.path.join(WApath,fname)
                        WA=np.genfromtxt(fdir, delimiter=",")
                        decode_precomp(well_assigner=WA, differentiate=diff, 
                        scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                        decname=os.path.join(dpath, 'decoder_'+method+'.json')
                        json.dump( dict_decode, open(decname, 'w' ) )
 


        N+=step




def decode_single( WA, readout:np.ndarray, 
                  differentiate, scrambler=True,
            dir_scramblers=False, **kwargs) -> list:
        print('Function not yet ready, come back later')
        return 'Function not yet ready, come back later'
        if scrambler==True:
            N=WA.shape[0]
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        elif not scrambler==False:
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        
        elif not dir_scramblers==False:
            N=WA.shape[0]
            if differentiate==1:
                scrambler={1:np.arange(N)}
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        

def str_to_tuple(string, delimiter='-'):
    return tuple(string.split(delimiter))

def tuple_to_str(tuple_type, delimiter='-'):
    return delimiter.join(map(str,tuple_type))
 
    



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



dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,
             'keep_ratios_constant': keep_ratios_constant,
             'start':start, 'stop':stop,  'step':step, 
             'dir_WAs':args.dir_WAs, 'dir_scramblers':args.dir_scramblers, 
             'max_differentiate': max_diff, 'max_dims':max_dims}


decode_sweep(**dict_kwargs)
