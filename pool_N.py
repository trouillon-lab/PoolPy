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
parser.add_argument('--guesses', type=int, default=5, help='An integer argument for guesses of random WA with default 20')
parser.add_argument('--method', type=str, default='all', help="A string argument with default 'all'")
parser.add_argument('--path', type=str, default='./pooling_results', help="A string argument with default './pooling_results'")

args = parser.parse_args()

args_dict = vars(args)

diff= args_dict['differentiate']
n_compounds=args_dict['n_compounds']
method=args_dict['method']
this_path=args_dict['path']
args_dict['save_dir']=copy.deepcopy(args_dict['path'])
args_dict['return_wa']=True
args_dict['keep_ratios_constant']=False


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
    WA_rand,  min_tests, perc_check=assign_wells_random_precomp(scrambler=scrambler, return_me=True, **args_dict )
    #thisfile=os.path.join(WA_path,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+
    #                            '_ME_'+str(np.round(min_tests,2))+'.csv')
    WA_list.append(WA_rand.astype(int))
    multi.append('Random')
    
if method=='matrix' or method=='all':
    WA_mat=assign_wells_mat(**args_dict)
    WA_list.append(WA_mat.astype(int))
    multi.append('Matrix')

if method=='binary' or method=='all':

    WA_bin=assign_wells_bin(**args_dict)
    multi.append('Binary')
    WA_list.append(WA_bin.astype(int))

if method=='std' or method=='all':

     WA_std=assign_wells_STD(**args_dict)
     multi.append('STD')
     WA_list.append(WA_std.astype(int))


if method=='chinese_trick' or method=='all':

    WA_chin=assign_wells_chinese(**args_dict)
    multi.append('Chinese trick')
    WA_list.append(WA_chin.astype(int))


this_dir=os.path.join(this_path,'N_'+str(args_dict['n_compounds']), 'diff_'+str(args_dict['differentiate']), 'WAs')

if not os.path.exists(this_dir):
    os.makedirs(this_dir)

for method, WA in zip(multi, WA_list):
    thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(args_dict['n_compounds'])+'_diff_'+str(args_dict['differentiate'])+'.csv')
    np.savetxt(thisfile, WA.astype(int), delimiter=",")
    


ls_names_met=['Method', 'Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps']
ls_met=[]
full_methods=[]
WApath=this_dir
filenames = next(os.walk(WApath), (None, None, []))[2]
for fname in filenames:
    #print(fname)
    fdir=os.path.join(WApath,fname)
    WA=np.genfromtxt(fdir, delimiter=",")
    mean_exp, extra_exp,  _, perc_check= mean_metrics_precomp(well_assigner=WA,scrambler=scrambler, **args_dict)
    n_wells=WA.shape[1]
    M_exp=np.round(mean_exp, 2)
    max_comp=np.max(np.sum(WA, axis=0))
    method=re.sub('^WA_', '', fname)
    method=re.sub('_.*$', '', method)
    #print(method)
    ls_met.append([method, M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100])
    full_methods.append(method)
Hier=calculate_metrics_hierarchical(**args_dict)
ls_met.append(['Hierarchical']+ [np.round(i,2) for i in Hier[:-1]])
full_methods.append('Hierarchical')
df_met=pd.DataFrame(ls_met)

this_dir=os.path.join(this_path,'N_'+str(args_dict['n_compounds']), 'diff_'+str(args_dict['differentiate']))

idx_renamer={i:j for i,j in zip(df_met.index, full_methods)}
col_renamer={i:j for i,j in zip(df_met.columns, ls_names_met)}
df_met.rename(index=idx_renamer, columns=col_renamer, inplace=True)
metname=os.path.join(this_dir, 'Metrics_N_'+str(n_compounds)+'_diff_'+str(diff)+'.csv')
df_met.to_csv(metname)

this_dir=os.path.join(this_path,'N_'+str(args_dict['n_compounds']), 'diff_'+str(args_dict['differentiate']), 'WAs')

for method, WA in zip(multi, WA_list):
    thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(args_dict['n_compounds'])+'_diff_'+str(args_dict['differentiate'])+'.csv')
    tmp1=pd.DataFrame(WA.astype(int), columns=['Pool '+ str(i) for i in range(WA.shape[1])], index=['Sample '+ str(i) for i in range(WA.shape[0])])
    tmp1.to_csv(thisfile)