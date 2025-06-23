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
import shutil



from Functions import *



def get_multidim_methods(n_compounds, max_dims):
    """Generate multidim method names based on current parameters"""
    multi = []
    for i in np.arange(2, int(np.ceil(np.log(n_compounds)/np.log(2)))):
        if i > max_dims:
            continue
        multi.append(f'multidim-{i}')
    return multi

def process_methods_for_kwargs(**kwargs):
    """
    Checks which WA files exist, computes missing ones, and copies diff-independent files
    to new diff folders as needed. Automatically includes multidim methods.
    """
    # Extract parameters with defaults
    n_compounds = kwargs['n_compounds']
    max_diff = kwargs['max_diff']
    save_dir = kwargs['save_dir']
    max_dims = kwargs.get('max_dims', 4)
    timeit = kwargs.get('timeit', False)
    
    # Generate method list dynamically
    multidim_methods = get_multidim_methods(n_compounds, max_dims)
    all_methods = multidim_methods + ['Matrix', 'Binary', 'STD', 'Chinese trick']
    
    # Categorize methods
    DIFF_INDEPENDENT = multidim_methods + ['Matrix', 'Binary']
    DIFF_DEPENDENT = ['STD', 'Chinese trick']

    for diff in range(1, max_diff + 1):
        diff_dir = os.path.join(save_dir, f'N_{n_compounds}', f'diff_{diff}', 'WAs')
        os.makedirs(diff_dir, exist_ok=True)  # Handles dir creation in one line

        for method in all_methods:
            # Construct filename using consistent pattern
            filename = f'WA_{method}_N_{n_compounds}_diff_{diff}.csv'
            wa_file = os.path.join(diff_dir, filename)
            
            if method in DIFF_INDEPENDENT:
                # For diff-independent methods
                if diff == 1:
                    if not os.path.exists(wa_file):
                        if timeit:
                            print(f"Computing {method} for diff=1")
                        # compute_and_save(method, n_compounds, diff, **kwargs)
                else:
                    # Check if diff=1 file exists to copy
                    diff1_file = os.path.join(
                        save_dir, 
                        f'N_{n_compounds}',
                        'diff_1',
                        'WAs',
                        filename.replace(f'_diff_{diff}', '_diff_1')
                    )
                    if not os.path.exists(wa_file) and os.path.exists(diff1_file):
                        shutil.copy(diff1_file, wa_file)
                        if timeit:
                            print(f"Copied {method} from diff=1 to diff={diff}")
            else:  # Diff-dependent methods
                if not os.path.exists(wa_file):
                    if timeit:
                        print(f"Computing {method} for diff={diff}")
                    # compute_and_save(method, n_compounds, diff, **kwargs)






def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    
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

    WA_bin=assign_wells_bin(**kwargs)
    methods.append('Binary')
    WA_list.append(WA_bin)

    if kwargs['max_diff']>1:
        for diffo in np.arange(kwargs['max_diff']):


            kwargs.update({'differentiate':int(diffo+1)})

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

            for method, WA in zip(methodso, WA_listo):
                thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_'+str(kwargs['differentiate'])+'.csv')
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



        this_dir=os.path.join(kwargs['save_dir'],'N_'+str(kwargs['n_compounds']), 'diff_1', 'WAs')

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        for method, WA in zip(methods, WA_list):
            thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_1.csv')
            np.savetxt(thisfile, WA.astype(bool), delimiter=",")









parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--save_dir')
parser.add_argument('--max_diff')
parser.add_argument('--max_dims')
parser.add_argument('--timeit')






args = parser.parse_args()



differentiate= 2 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
step= 10 if type(args.step)==type(None) else int(args.step)
save_dir= os.path.join(os.getcwd(),'outs') if type(args.save_dir)==type(None) else str(args.save_dir)
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'
max_diff= 4 if type(args.max_diff)==type(None) else int(args.max_diff)
max_dims= 4 if type(args.max_dims)==type(None) else int(args.max_dims)



dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,
             'start':start, 'stop':stop,  'step':step, 'save_dir':save_dir, 'max_diff': max_diff, 'max_dims':max_dims}


make_all_deterministic_WAs(**dict_kwargs)