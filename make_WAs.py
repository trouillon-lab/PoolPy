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

import numpy as np
import os
import shutil
import copy
import time
from Functions import assign_wells_multidim, assign_wells_mat, assign_wells_bin, assign_wells_STD, assign_wells_chinese

def get_wa_filename(save_dir, n_compounds, diff, method):
    """Generate consistent filename pattern"""
    return os.path.join(
        save_dir,
        f'N_{n_compounds}',
        f'diff_{diff}',
        'WAs',
        f'WA_{method}_N_{n_compounds}_diff_{diff}.csv'
    )

def process_n_compounds(**kwargs):
    """
    Processes WA computations for a specific n_compounds value with diff handling
    and file existence checks
    """
    n_compounds = kwargs['n_compounds']
    max_diff = kwargs['max_diff']
    save_dir = kwargs['save_dir']
    max_dims = kwargs['max_dims']
    timeit = kwargs.get('timeit', False)
    
    # Generate multidim methods dynamically
    multidim_methods = []
    for i in np.arange(2, int(np.ceil(np.log(n_compounds)/np.log(2)))):
        if i > max_dims:
            continue
        multidim_methods.append(f'multidim-{i}')
    
    # Compute diff-independent methods only if needed
    WA_list = []
    methods = []
    computed_diff_independent = False
    
    # Check if all diff-independent files exist for diff=1
    all_exist = True
    for method in multidim_methods + ['Matrix', 'Binary']:
        if not os.path.exists(get_wa_filename(save_dir, n_compounds, 1, method)):
            all_exist = False
            break
    
    # Compute diff-independent methods if any are missing
    if not all_exist:
        computed_diff_independent = True
        # 1. Compute multidim methods
        for method in multidim_methods:
            dim = int(method.split('-')[1])
            WA = assign_wells_multidim(n_dims=dim, **kwargs)
            WA_list.append(WA)
            methods.append(method)
        
        # 2. Compute Matrix and Binary
        WA_mat = assign_wells_mat(**kwargs)
        WA_list.append(WA_mat)
        methods.append('Matrix')
        
        WA_bin = assign_wells_bin(**kwargs)
        WA_list.append(WA_bin)
        methods.append('Binary')
    elif timeit:
        print(f"All diff-independent files exist for n={n_compounds}, skipping computation")
    
    # Process each differentiation level
    for diff in range(1, max_diff + 1):
        current_kwargs = kwargs.copy()
        current_kwargs['differentiate'] = diff
        
        # Create directory structure
        diff_dir = os.path.join(save_dir, f'N_{n_compounds}', f'diff_{diff}', 'WAs')
        os.makedirs(diff_dir, exist_ok=True)
        
        # Handle diff-independent methods
        for method in multidim_methods + ['Matrix', 'Binary']:
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                if diff == 1 and computed_diff_independent:
                    # Find method index and save
                    idx = methods.index(method)
                    np.savetxt(dst_file, WA_list[idx].astype(bool), delimiter=",")
                    if timeit:
                        print(f"Saved {method} for diff={diff}")
                else:
                    # Try to copy from diff=1
                    src_file = get_wa_filename(save_dir, n_compounds, 1, method)
                    if os.path.exists(src_file):
                        shutil.copy(src_file, dst_file)
                        if timeit:
                            print(f"Copied {method} to diff={diff}")
                    elif timeit:
                        print(f"Warning: Base file missing for {method} at diff={diff}")
        
        # Handle diff-dependent methods
        for method in ['STD', 'Chinese trick']:
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                # Compute only if file doesn't exist
                if method == 'STD':
                    WA = assign_wells_STD(**current_kwargs)
                else:  # Chinese trick
                    WA = assign_wells_chinese(**current_kwargs)
                
                np.savetxt(dst_file, WA.astype(bool), delimiter=",")
                if timeit:
                    print(f"Computed {method} for diff={diff}")
            elif timeit:
                print(f"Skipping {method} for diff={diff} (already exists)")

def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    """
    Main loop to process all n_compounds values
    """
    current = start
    while current < stop:
        if kwargs.get('timeit'):
            print(f"Processing {current} compounds")
            time0 = time.time()
        
        # Process this n_compounds value
        kwargs['n_compounds'] = current
        process_n_compounds(**kwargs)
        
        if kwargs.get('timeit'):
            DTS=np.round((time.time() - time0),2)
            DTD=DTS//86400
            DTH=DTS//3600-DTD*24
            DTM=DTS//60-DTH*60-DTD*24*60
            DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
            print('\n')
            print("%s days %s hours %s minutes and %s seconds required for N= %s compounds" % 
                  (DTD, DTH, DTM, DTS, current))
            print('\n')

        
        current += step

# Argument parsing and main execution remains the same as in your original code
# ... [rest of your argument parsing and main call] ...




def make_all_deterministic_WAs_old(start=50, stop=150, step=10, **kwargs):
    
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