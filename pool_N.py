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
from Fast_functions import *
import argparse

parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=2, help='An integer argument with default 2')
parser.add_argument('--n_compounds', type=int, default=50, help='An integer argument with default 50')
parser.add_argument('--method', type=str, default='all', help="A string argument with default 'all'")
parser.add_argument('--path', type=str, default='./', help="A string argument with default './'")

args = parser.parse_args()

args_dict = vars(args)

diff= args_dict['differentiate']
n_compounds=args_dict['n_compounds']
method=args_dict['method']
this_path=args_dict['path']
args_dict['save_dir']=copy.deepcopy(args_dict['path'])
WA_path=os.path.join(this_path, 'WAs')

scrambler={1:np.arange(n_compounds)}
for j in range(2,diff+1):
    scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
WA_list=[]
multi=[]

if method.startswith('multidim') or method=='all':

    for i in np.arange(2,int(np.ceil(np.log(n_compounds)/np.log(2)))):
        WA_mul=assign_wells_multidim(n_dims=i, **args_dict)
        WA_list.append(WA_mul)
        multi.append('multidim-'+str(i))

if method=='random' or method=='all':
    WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                                differentiate=diff,scrambler=scrambler, return_me=True, **args_dict )
    thisfile=os.path.join(WA_path,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+
                                '_ME_'+str(np.round(min_tests,2))+'.csv')
    
if method=='matrix' or method=='all':
    WA_mat=assign_wells_mat(**args_dict)
    WA_list.append(WA_mat)
    multi.append('Matrix')

if method=='binary' or method=='all':

    WA_bin=assign_wells_bin(**args_dict)
    multi.append('Binary')
    WA_list.append(WA_bin)

if method=='std' or method=='all':

     WA_std=assign_wells_STD(**args_dict)
     multi.append('STD')
     WA_list.append(WA_std)


if method=='chinese_trick' or method=='all':

    WA_chin=assign_wells_chinese(**args_dict)
    multi.append(''Chinese trick'')
    WA_list.append(WA_chin)