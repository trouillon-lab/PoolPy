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
parser.add_argument('--differentiate', type=int, default=-1, help='An integer argument with default 2')
parser.add_argument('--path_to_WA', type=str, help="A string argument with default './pooling_results'")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout")

args = parser.parse_args()

args_dict = vars(args)



dira=args_dict['path_to_WA']
diff=args_dict['differentiate']
readout_in=args_dict['readout']
WA_df=pd.read_csv(dira, index_col=0)

if readout_in.endswith('csv'):
    readout = np.genfromtxt('readout_in', delimiter=',', dtype=int)
else:
    readout = np.fromstring(readout_in, sep=',', dtype=int)    

st_dir=str(dira)
inf_diff=re.sub('^.*_', '', st_dir)
inf_diff=re.sub('\\.csv$', '', inf_diff)
inf_diff=int(inf_diff)

method=re.sub('.*WA_', '', dira)
method=re.sub('_.*$', '', method)

if diff==-1:
    diff=inf_diff

if diff!=inf_diff:
    print(f'WARNING: inferred differentiate of {inf_diff} different from passed differentiate of {diff}\n')


WA=WA_df.values
n_compounds=WA.shape[1]
scrambler={1:np.arange(n_compounds)}
for j in range(2,diff+1):
    scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
WA_list=[]
multi=[]

scrambler={1:np.arange(n_compounds)}
for j in range(2,diff+1):
    scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
WA_list=[]
multi=[]

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

    f.write('Decoded file {dira} assuming differentiate {diff}.\n Possible positive samples combiantion are')
    for line in decoded:
        f.write(f"Samples: {line}\n")
