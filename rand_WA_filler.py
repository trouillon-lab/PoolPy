import numpy as np
import os
import shutil
import time
import argparse
from Functions import assign_wells_random_precomp  # Import your random method function

def process_random_method(n_compounds, diff, save_dir, **kwargs):
    """Process random method for a specific n_compounds and diff"""
    timeit = kwargs.get('timeit', False)
    WApath = os.path.join(save_dir, f'N_{n_compounds}', f'diff_{diff}', 'WAs')
    os.makedirs(WApath, exist_ok=True)
    
    # Build filename with metrics (same pattern as your code)
    filename = f'WA_Random_N_{n_compounds}_diff_{diff}_ME_0.0.csv'  # Placeholder for metrics
    
    # Check if file already exists
    if os.path.exists(os.path.join(WApath, filename)) and not kwargs.get('overwrite', False):
        if timeit:
            print(f"Skipping Random method for N={n_compounds}, diff={diff} (file exists)")
        return
    
    # Compute random well assignment
    scrambler = {1: np.arange(n_compounds)}
    WA_rand, min_tests, perc_check = assign_wells_random_precomp(
        n_compounds=n_compounds,
        differentiate=diff,
        scrambler=scrambler,
        return_me=True,
        **kwargs
    )
    
    # Update filename with actual metrics
    filename = f'WA_Random_N_{n_compounds}_diff_{diff}_ME_{np.round(min_tests,2)}.csv'
    output_file = os.path.join(WApath, filename)
    
    # Save results
    np.savetxt(output_file, WA_rand.astype(bool), delimiter=",")
    if timeit:
        print(f"Computed Random method for N={n_compounds}, diff={diff}")

def process_n_compounds(**kwargs):
    """Process all methods for a specific n_compounds value"""
    n_compounds = kwargs['n_compounds']
    max_diff = kwargs['max_diff']
    save_dir = kwargs['save_dir']
    timeit = kwargs.get('timeit', False)
    
    # Process each differentiation level
    for diff in range(1, max_diff + 1):
        if timeit:
            print(f"Processing N={n_compounds}, diff={diff}")
        
        # Process deterministic methods (your existing code)
        # ... [your deterministic processing code] ...
        
        # Process random method
        process_random_method(n_compounds, diff, save_dir, **kwargs)

def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    """Main loop to process all n_compounds values"""
    current = start
    while current < stop:
        if kwargs.get('timeit'):
            print(f"Processing n_compounds={current}")
            time0 = time.time()
        
        # Process this n_compounds value
        kwargs['n_compounds'] = current
        process_n_compounds(**kwargs)
        
        if kwargs.get('timeit'):
            elapsed = np.round(time.time() - time0, 1)
            print(f"Completed n_compounds={current} in {elapsed} seconds")
        
        current += step

# Argument parsing remains similar to your original code
parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--save_dir')
parser.add_argument('--max_diff')
parser.add_argument('--max_dims')
parser.add_argument('--timeit')
parser.add_argument('--overwrite')

args = parser.parse_args()

# Set parameters with defaults
differentiate = 2 if args.differentiate is None else int(args.differentiate)
start = 50 if args.start is None else int(args.start)
stop = 110 if args.stop is None else int(args.stop)
step = 10 if args.step is None else int(args.step)
save_dir = os.path.join(os.getcwd(), 'outs') if args.save_dir is None else args.save_dir
timeit = True if args.timeit is None else args.timeit == 'True'
max_diff = 4 if args.max_diff is None else int(args.max_diff)
max_dims = 4 if args.max_dims is None else int(args.max_dims)
overwrite = False if args.overwrite is None else args.overwrite == 'True'

dict_kwargs = {
    'differentiate': differentiate,
    'return_wa': True,
    'timeit': timeit,
    'start': start,
    'stop': stop,
    'step': step,
    'save_dir': save_dir,
    'max_diff': max_diff,
    'max_dims': max_dims,
    'overwrite': overwrite
}

make_all_deterministic_WAs(**dict_kwargs)
