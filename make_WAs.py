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



def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    dict_comp={}
    current=start
    while current<stop:
        time0=time.time()
        if kwargs['timeit']:
            print(current)
        full_deterministic_WAS(n_compounds=current, **kwargs)
        current=current+step
        if kwargs['timeit']:
            print("segment time: %s seconds" % np.round(time.time() - time0, 1))


def full_deterministic_WAS(**kwargs):
    #methods=['matrix', 'random', 'STD', 'Chinese trick']
    # matrix assignment
    
    kwargs['return_wa']=True

    WA_list=[]
    multi=[]
    for i in np.arange(2,int(np.ceil(np.log(kwargs['n_compounds'])/np.log(2)))):
        if i>kwargs['max_dims']:
            continue
        WA_mul=assign_wells_multidim(n_dims=i, **kwargs)
        WA_list.append(WA_mul)
        multi.append('multidim-'+str(i))
        
    methods=multi.copy()

    WA_mat=assign_wells_mat(**kwargs)

    WA_list.append(WA_mat)
    methods.append('Matrix')

    WA_bin=assign_wells_L(**kwargs)
    methods.append('Binary')
    WA_list.append(WA_bin)

    if kwargs['max_diff']>1:
        for diffo in np.arange(kwargs['max_diff']):
            diff=diffo+1

            WA_listo=copy.deepcopy(WA_list)
            methodso=copy.deepcopy(methods)

            # STD asignment 
            WA_std=assign_wells_STD(differentiate=diff,**kwargs)

            

            # chinese trick assignment
            WA_chin=assign_wells_chinese(differentiate=diff,**kwargs)


            WA_listo.extend([ WA_std, WA_chin])
            methodso.extend(['STD', 'Chinese trick'])



            this_dir=os.path.join(kwargs['save_dir'],'N_'+str(kwargs['n_compounds']), 'diff_'+str(diff), 'WAs')

            if not os.path.exists(this_dir):
                os.makedirs(this_dir)

            for method, WA in zip(methods, WA_list):
                thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_'+str(diff)+'.csv')
                np.savetxt(thisfile, WA.astype(bool), delimiter=",")
    
    else:
        


        WA_listo=copy.deepcopy(WA_list)
        methodso=copy.deepcopy(methods)

        # STD asignment 
        WA_std=assign_wells_STD(**kwargs)

        

        # chinese trick assignment
        WA_chin=assign_wells_chinese(**kwargs)


        WA_listo.extend([ WA_std, WA_chin])
        methodso.extend(['STD', 'Chinese trick'])



        this_dir=os.path.join(kwargs['save_dir'],'N_'+str(kwargs['n_compounds']), 'diff_'+str(kwargs['differentiate']), 'WAs')

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        for method, WA in zip(methods, WA_list):
            thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_'+str(kwargs['differentiate']+'.csv'))
            np.savetxt(thisfile, WA.astype(bool), delimiter=",")









parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--base_dir')
parser.add_argument('--max_diff')
parser.add_argument('--timeit')






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