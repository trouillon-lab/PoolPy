import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import pickle
import time
import os


from Functions import *


parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--base_dir')
parser.add_argument('--rand_guesses')
parser.add_argument('--max_compounds')
parser.add_argument('--n_compounds_per_well')
parser.add_argument('--n_wells')
parser.add_argument('--n_dims')
parser.add_argument('--return_wa')
parser.add_argument('--timeit')
parser.add_argument('--method')
parser.add_argument('--inline_print')
parser.add_argument('--keep_ratios_constant')
parser.add_argument('--all_dims')
parser.add_argument('--max_dims')





args = parser.parse_args()



differentiate= 2 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
step= 10 if type(args.step)==type(None) else int(args.step)
base_dir= os.getcwd() if type(args.base_dir)==type(None) else str(args.base_dir)
rand_guesses= 10 if type(args.rand_guesses)==type(None) else int(args.rand_guesses)
return_wa= True if type(args.return_wa)==type(None) else args.return_wa=='True'
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'
inline_print= False if type(args.inline_print)==type(None) else args.inline_print=='True'
keep_ratios_constant= False if type(args.keep_ratios_constant)==type(None) else args.keep_ratios_constant=='True' 
all_dims= False if type(args.all_dims)==type(None) else args.all_dims=='True'
max_dims= np.inf if type(args.max_dims)==type(None) else int(args.max_dims)

dict_kwargs={'differentiate':differentiate, 'return_wa':return_wa, 'timeit':timeit,
             'start':start, 'stop':stop,  'step':step, 'base_dir':base_dir, 'rand_guesses':rand_guesses,
             'inline_print':inline_print, 'keep_ratios_constant': keep_ratios_constant, 'all_dims':all_dims, 'max_dims':max_dims}

if type(args.max_compounds)!=type(None): 
    dict_kwargs.update({'max_compounds':int(args.max_compounds)})
if type(args.n_compounds_per_well)!=type(None): 
    dict_kwargs.update({'n_compounds_per_well':int(args.n_compounds_per_well)})
if type(args.n_wells)!=type(None): 
    dict_kwargs.update({'n_wells':int(args.n_wells)})
if type(args.n_dims)!=type(None): 
    dict_kwargs.update({'n_dims':int(args.n_dims)})
if type(args.method)!=type(None): 
    dict_kwargs.update({'method':args.method})

print(dict_kwargs)

start_time = time.time()

if 'method' in dict_kwargs.keys():
    if timeit:
        print('\n')
        print('-----------------------------------------------------')
        print("total time: %s seconds" % np.round(time.time() - start_time, 1))
        print('-----------------------------------------------------')
    dict_c=single_method_sweep(**dict_kwargs)
    if not dict_kwargs['inline_print'] or dict_kwargs['method'] not in ['STD', 'CT']:
        dict_c.update({'kwargs':dict_kwargs})


        
        method=dict_kwargs['method']
        if 'n_dims' in dict_kwargs.keys() and dict_kwargs['method']=='multidim':
            method=method+'_'+str(dict_kwargs['n_dims'])


        fpath=os.path.join(base_dir,'single_method',method)

        if not os.path.exists(fpath):
            os.makedirs(fpath)

        full_dir=os.path.join(fpath,method+'__'+str(start)+'-'+str(stop)+'_step'+str(step)+'.pk')
        with open(full_dir, 'wb') as handle:
            pickle.dump(dict_c, handle, protocol=pickle.HIGHEST_PROTOCOL)




else:
    dict_c=full_sweep_comparison(**dict_kwargs)
    if timeit:
        print('\n')
        print('-----------------------------------------------------')
        print("total time: %s seconds" % np.round(time.time() - start_time, 1))
        print('-----------------------------------------------------')



    dict_c.update({'kwargs':dict_kwargs})

    fpath=os.path.join(base_dir,'diff_'+str(differentiate))

    if not os.path.exists(fpath):
        os.makedirs(fpath)

    full_dir=os.path.join(fpath,str(start)+'-'+str(stop)+'_step'+str(step)+'.pk')



    with open(full_dir, 'wb') as handle:
        pickle.dump(dict_c, handle, protocol=pickle.HIGHEST_PROTOCOL)



# How to use :
# python pre-computation.py --start 20 --stop 47 --step 5 --differentiate 1 --rand_guesses 3 --base_dir 'where/you/want' 
# other inputs which are relevant only for random are --max_compounds --n_compounds_per_well --n_wells
# only for multidim --n_dims