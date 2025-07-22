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
parser.add_argument('--differentiate', type=int, default=2, help='An integer argument with default 2')
parser.add_argument('--path_to_WA', type=str, help="A string argument with default './pooling_results'")

args = parser.parse_args()

args_dict = vars(args)

dira=args_dict['path_to_WA']
diff=args_dict['differentiate']
WA_df=pd.read_csv(dira, index_col=0)

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
               readout=np.array((WA[0]+WA[5]).astype(bool).astype(int)))

print('The possible positives for the given well assigner, outcome, and differentiate are:')
for deco in decoded:
    print('Samples:', deco)


with open('your_file.txt', 'w') as f:
    for line in decoded:
        f.write(f"Samples: {line}\n")
