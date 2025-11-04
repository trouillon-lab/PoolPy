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
from Fast_functions import *



parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=-1, help='The number of maximum number of positives in the pooling strategy.')
parser.add_argument('--path_to_WA', type=str, help="A string argument containing the path to the pooling strategy file.")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout (readout in the form 0,1,0,1,0,0 or 3,6,14)")
parser.add_argument('--extensive_search', type=str, default='False', help="weather to search all possibilities (True) or use a faster exclusion based implementation (False, default)")

args = parser.parse_args()

args_dict = vars(args)



dira=args_dict['path_to_WA']
diff=args_dict['differentiate']
readout_in=args_dict['readout']
search_all=args_dict['extensive_search'].lower()=='true'
WA_df=pd.read_csv(dira, index_col=0)

if readout_in.endswith('csv'):
    readout = pd.read_csv(readout_in,header=None).to_numpy().reshape(-1)
else:
    readout = np.fromstring(readout_in, sep=',', dtype=int)   

WA=WA_df.values
n_pools=WA.shape[1]
n_compounds=WA.shape[0]

if np.max(readout)>1 or len(readout)!=n_pools:
    readout_bin_ls = [1 if i in readout else 0 for i in range(n_pools)]
    readout=np.array(readout_bin_ls)
 

st_dir=str(dira)
#inf_diff=re.sub('^.*_', '', st_dir)
pth_wa=re.sub('\\.csv$', '', st_dir)
#inf_diff=int(inf_diff)

#method=re.sub('.*WA_', '', dira)
#method=re.sub('_.*$', '', method)
readout_bl=np.array(readout.astype(bool).astype(int))

if search_all:
    if diff==-1:
        diff=n_compounds

    scrambler={1:np.arange(n_compounds)}
    for j in range(2,diff+1):
        scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
        

    decoded=decode_precomp(well_assigner=WA, differentiate= diff, scrambler=scrambler, 
                readout=readout_bl)
        

    if len(decoded)==0:
        print('We found no matches for the given parameters, check your input or try increasing the differentiate value')


    elif len(decoded)>n_compounds:
        print('The possible combinations of positive samples resulting in your readout is too large. Either test them individually or change pooling strategy')

    else:
        print('The possible positives for the given well assigner, outcome, and differentiate are:')
        for deco in decoded:
            print('Samples:', deco)
        
    

else: 
    mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
    original_indices = np.where(mask)[0]  # Get original row indices
    filtered_WA = WA[mask]
    n_compounds=filtered_WA.shape[0]
    if diff==-1:
        diff=n_compounds

    if n_compounds<2:
        if n_compounds==1:
            decoded=[original_indices[0]]
            print('The possible positives for the given well assigner, outcome, and differentiate are:')
            for deco in decoded:
                print('Samples:', deco)
        else:
            print('We found no matches for the given parameters, check your input or try increasing the differentiate value')


    else:

        scrambler={1:np.arange(n_compounds)}
        for j in range(2,diff+1):
            scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
        decoded_pre=decode_precomp(well_assigner=filtered_WA, differentiate= diff, scrambler=scrambler, 
                readout=readout_bl)
        
        # Map filtered indices back to original indices
        decoded = [[int(original_indices[idx]) for idx in combination] for combination in decoded_pre]
        
        # Remove duplicate combinations (convert to tuples for uniqueness, then back to lists)

        if len(decoded)==0:
            print('We found no matches for the given parameters, check your input or try increasing the differentiate value')


        elif len(decoded)>n_compounds:
            decoded_set = list(set([x for combo in decoded for x in combo]))
            decoded_set.sort()
            print('The possible combinations of positive samples resulting in your readout is too large.\n Either test all putative positives individually or change pooling strategy')
            print('The set of putative positives is ', decoded_set)

        else:
            print('The possible positives for the given well assigner, outcome, and differentiate are:')
            for deco in decoded:
                print('Samples:', deco)









fdriro=pth_wa+'_decoded.txt'

with open(fdriro, 'w+') as f:
    if len(decoded)>n_compounds:
        if search_all:
            f.write('The possible combinations of positive samples resulting in your readout is too large.\n Either test them individually or change pooling strategy')
        else:
            f.write('The possible combinations of positive samples resulting in your readout is too large.\n Either test all putative positives individually or change pooling strategy\n')
            f.write(f'The set of putative positives is {decoded_set}\n')
        exit()

    f.write(f'Decoded reaout {[i for i in range(n_pools) if readout[i]==1]} file {dira} with maximumn {diff} positives.\n Possible positive samples combinations are\n')
    for line in decoded:
        f.write(f"Samples: {line}\n")
