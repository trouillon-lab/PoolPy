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
parser.add_argument('--differentiate', type=int, default=-1, help='An integer argument for the differentiate value. Is inferred from WA file name by default.')
parser.add_argument('--path_to_WA', type=str, help="A string argument containing the path to the well assigner. WA file name should follow format given by app or pool_N.py.")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout (readout in the form 0,1,0,1,0,0 or 3,6,14)")

args = parser.parse_args()

args_dict = vars(args)



dira=args_dict['path_to_WA']
diff=args_dict['differentiate']
readout_in=args_dict['readout']
WA_df=pd.read_csv(dira, index_col=0)

if readout_in.endswith('csv'):
    readout = pd.read_csv(readout_in,header=None).to_numpy().reshape(-1)
else:
    readout = np.fromstring(readout_in, sep=',', dtype=int)   

WA=WA_df.values
n_pools=WA.shape[1]

if np.max(readout)>1 or len(readout)!=n_pools:
    readout_bin_ls = [1 if i in readout else 0 for i in range(n_pools)]
    readout=np.array(readout_bin_ls)
 

st_dir=str(dira)
inf_diff=re.sub('^.*_', '', st_dir)
inf_diff=re.sub('\\.csv$', '', inf_diff)
inf_diff=int(inf_diff)

method=re.sub('.*WA_', '', dira)
method=re.sub('_.*$', '', method)

if diff==-1:
    diff=inf_diff

if diff!=inf_diff:
    print(f'WARNING: differentiate of {inf_diff} inferred from WA file name different from passed differentiate of {diff}\n')

scrambler={1:np.arange(n_pools)}
for j in range(2,diff+1):
    scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_pools),j)))})
    

decoded=decode_precomp(well_assigner=WA,differentiate= diff, scrambler=scrambler, 
               readout=np.array(readout.astype(bool).astype(int)))

print('The possible positives for the given well assigner, outcome, and differentiate are:')
for deco in decoded:
    print('Samples:', deco)

s1driro=os.path.join(os.path.dirname(os.path.dirname(dira)), 'decoded')
if not os.path.isdir(s1driro):
    os.mkdir(s1driro)
fdriro=os.path.join(s1driro, f'{method}_diff_{diff}_decoded.txt')

with open(fdriro, 'w+') as f:
    if diff!=inf_diff:
        f.write(f'WARNING: inferred differentiate of {inf_diff} different from passed differentiate of {diff}\n')

    f.write(f'Decoded file {dira} assuming differentiate {diff}.\n Possible positive samples combinations are ')
    for line in decoded:
        f.write(f"Samples: {line}\n")
