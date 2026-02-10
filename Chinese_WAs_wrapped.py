import numpy as np
import os
import time
import math
import itertools
import string
from Functions_wrapped import (
    assign_wells_chinese,
    format_time,
    int_to_base
)


def get_wa_filename(save_dir, n_compounds, diff, method):
    """Generate consistent filename pattern"""
    return os.path.join(
        save_dir,
        f'N_{n_compounds}',
        f'diff_{diff}',
        'WAs',
        f'WA_{method}_N_{n_compounds}_diff_{diff}.csv'
    )


def remove_chinese_wa_file(filepath):
    """Remove a WA file if it exists"""
    if os.path.exists(filepath):
        try:
            os.remove(filepath)
            return True
        except Exception as e:
            print(f"Warning: Could not remove {filepath}: {e}")
            return False
    return False


def process_chinese_wa(n_compounds, diff, save_dir, method_type, timeit=False, remove_existing=True):
    """
    Process and save Chinese remainder WA for given parameters
    
    Parameters:
    -----------
    n_compounds : int
        Number of compounds
    diff : int
        Differentiation level
    save_dir : str
        Base directory where files will be saved
    method_type : str
        One of 'normal', 'special', 'backtrack'
    timeit : bool
        Whether to print timing information
    remove_existing : bool
        If True, remove the specific WA file being regenerated before computing
    """
    
    # Create directory structure
    diff_dir = os.path.join(save_dir, f'N_{n_compounds}', f'diff_{diff}', 'WAs')
    os.makedirs(diff_dir, exist_ok=True)
    
    # Determine method name and kwargs
    if method_type == 'normal':
        method_name = 'Chinese remainder'
        kwargs = {'special_diff': False, 'backtrack': False}
    elif method_type == 'special':
        method_name = 'Chinese special'
        kwargs = {'special_diff': True, 'backtrack': False}
    elif method_type == 'backtrack':
        method_name = 'Ch. rm. bktrk'
        kwargs = {'special_diff': False, 'backtrack': True}
    else:
        raise ValueError(f"Invalid method_type: {method_type}. Must be 'normal', 'special', or 'backtrack'")
    
    # Generate filename
    dst_file = get_wa_filename(save_dir, n_compounds, diff, method_name)
    
    # Remove existing file only for this specific (n_compounds, diff, method) combination
    if remove_existing:
        removed = remove_chinese_wa_file(dst_file)
        if timeit and removed:
            print(f"Removed existing {method_name} WA file for N={n_compounds}, diff={diff}")
    
    # Compute WA
    try:
        WA = assign_wells_chinese(
            n_compounds=n_compounds,
            differentiate=diff,
            **kwargs
        )
        
        # Save to file
        np.savetxt(dst_file, WA.astype(int), delimiter=",", fmt='%d')
        
        if timeit:
            print(f"Computed and saved {method_name} for N={n_compounds}, diff={diff}")
        
        return True
    except Exception as e:
        print(f"Error computing {method_name} for N={n_compounds}, diff={diff}: {e}")
        return False


def make_chinese_WAs(start=50, stop=150, step=10, 
                     start_diff=1, end_diff=None, 
                     method_type='all', save_dir='./wrapped',
                     max_prev=0.1, timeit=True, remove_existing=True, **kwargs):
    """
    Generate Chinese remainder WAs for a range of n_compounds and differentiation values.
    
    Parameters:
    -----------
    start : int
        Starting value for n_compounds
    stop : int
        Stopping value for n_compounds (inclusive)
    step : int
        Step size for n_compounds increments
    start_diff : int
        Starting value for differentiation
    end_diff : int or None
        Ending value for differentiation (inclusive). If None, computed as int(max_prev*n_compounds)+1
    method_type : str
        One of 'normal', 'special', 'backtrack', or 'all'
    save_dir : str
        Base directory where WA files will be saved
    max_prev : float
        Maximum differentiation as fraction of n_compounds
    timeit : bool
        Whether to print timing information
    remove_existing : bool
        If True, remove only the specific WA files being regenerated (default True).
        Only data for the current (n_compounds, diff, method) will be removed, not other Chinese WA files.
    """
    
    # Validate method_type
    valid_types = ['normal', 'special', 'backtrack', 'all']
    if method_type not in valid_types:
        raise ValueError(f"Invalid method_type: {method_type}. Must be one of {valid_types}")
    
    # Determine which methods to compute
    if method_type == 'all':
        methods_to_compute = ['normal', 'special', 'backtrack']
    else:
        methods_to_compute = [method_type]
    
    # Overall timer
    overall_start = time.time()
    
    # Loop over n_compounds values
    n_compounds = start
    while n_compounds <= stop:
        n_start_time = time.time()
        
        if timeit:
            print(f"\n{'='*60}")
            print(f"Processing N = {n_compounds} compounds")
            print(f"{'='*60}")
        
        # Determine max differentiation for this n_compounds
        max_diff = end_diff if end_diff is not None else int(max_prev * n_compounds) + 1
        
        # Loop over differentiation values
        for diff in range(start_diff, max_diff + 1):
            if diff > max_prev * n_compounds:
                break
            
            # Process each method type
            for method in methods_to_compute:
                # Skip 'special' method if diff not in [2, 3]
                if method == 'special' and diff not in [2, 3]:
                    if timeit:
                        print(f"Skipping 'special' method for diff={diff} (only valid for diff=2 or 3)")
                    continue
                
                process_chinese_wa(
                    n_compounds=n_compounds,
                    diff=diff,
                    save_dir=save_dir,
                    method_type=method,
                    timeit=timeit,
                    remove_existing=remove_existing
                )
        
        # Print time elapsed for this N
        n_elapsed = time.time() - n_start_time
        if timeit:
            print(f"\nCompleted N={n_compounds}, time elapsed {format_time(n_elapsed)}")
            print(f"{'='*60}\n")
        
        n_compounds += step
    
    # Print total completion time
    total_elapsed = time.time() - overall_start
    if timeit:
        print(f"\n{'#'*60}")
        print(f"All Chinese remainder WAs completed!")
        print(f"Total time elapsed: {format_time(total_elapsed)}")
        print(f"{'#'*60}\n")


if __name__ == '__main__':
    # Example usage
    make_chinese_WAs(
        start=10,
        stop=1000,
        step=1,
        start_diff=2,
        end_diff=3,
        method_type='special',
        save_dir='./wrapped',
        max_prev=0.1,
        timeit=True
    )
