import numpy as np
import argparse
import time
import os

from Fast_functions import *


parser = argparse.ArgumentParser()
parser.add_argument('--differentiate', type=int, default=2)
parser.add_argument('--start', type=int, default=50)
parser.add_argument('--stop', type=int, default=110)
parser.add_argument('--step', type=int, default=10)
parser.add_argument('--directory', type=str, default='./wrapped')
parser.add_argument('--max_diff', type=int, default=10)
parser.add_argument('--max_prev', type=float, default=0.1)
parser.add_argument('--max_dims', type=int, default=np.inf)
parser.add_argument('--timeit', type=str, default='True')
parser.add_argument('--cleanup', type=str, default='False')
parser.add_argument('--one_liner', type=str, default='True')
parser.add_argument('--max_compounds', type=int, default=None)
parser.add_argument('--n_compounds_per_well', type=int, default=None)
parser.add_argument('--n_wells', type=int, default=None)
parser.add_argument('--rand_guesses', type=int, default=10)
parser.add_argument('--max_redundancy', type=float, default=2)
parser.add_argument('--min_redundancy', type=float, default=0.5)
parser.add_argument('--overwrite', type=str, default='True')
parser.add_argument('--n_samp', type=int, default=0)

args = parser.parse_args()

start_time = time.time()

if  args.n_samp>0:
    dict_kwargs = {
        'differentiate': args.differentiate,
        'return_wa': True,
        'timeit': args.timeit == 'True',
        'guesses': args.rand_guesses,
        'start': args.n_samp,
        'stop': args.n_samp+1,
        'step': 2,
        'dir_WAs': args.directory,
        'save_dir': args.directory,
        'max_diff': args.max_diff,
        'max_dims': args.max_dims,
        'max_redundancy': args.max_redundancy,
        'min_redundancy': args.min_redundancy,
        'one_liner': args.one_liner == 'True',
        'overwrite': args.overwrite == 'True',
        'max_prev':args.max_prev,
        'cleanup': args.cleanup == 'True'
    }

else:
    dict_kwargs = {
        'differentiate': args.differentiate,
        'return_wa': True,
        'timeit': args.timeit == 'True',
        'guesses': args.rand_guesses,
        'start': args.start,
        'stop': args.stop,
        'step': args.step,
        'dir_WAs': args.directory,
        'save_dir': args.directory,
        'max_diff': args.max_diff,
        'max_dims': args.max_dims,
        'max_redundancy': args.max_redundancy,
        'min_redundancy': args.min_redundancy,
        'one_liner': args.one_liner == 'True',
        'overwrite': args.overwrite == 'True',
        'max_prev':args.max_prev,
        'cleanup': args.cleanup == 'True'
    }


make_all_deterministic_WAs(**dict_kwargs)

DTS = np.round((time.time() - start_time), 2)
DTD = DTS // 86400
DTH = DTS // 3600 - DTD * 24
DTM = DTS // 60 - DTH * 60 - DTD * 24 * 60
DTS = np.round(DTS - (DTM + DTH * 60 + DTD * 24 * 60) * 60, 2)

print('\n')
print('\n')
print('**********************************************************************************************************')
print('**********************************************************************************************************')
print('**********************************************************************************************************')
print('\n')
print("%s days %s hours %s minutes and %s seconds overall required for deterministic WAs" % (DTD, DTH, DTM, DTS))
print('\n')
print('\n')


# Add optional arguments if provided
if args.max_compounds is not None:
    dict_kwargs['max_compounds'] = args.max_compounds
if args.n_compounds_per_well is not None:
    dict_kwargs['n_compounds_per_well'] = args.n_compounds_per_well
if args.n_wells is not None:
    dict_kwargs['n_wells'] = args.n_wells

rand_N_sweep(**dict_kwargs)

DTS = np.round((time.time() - start_time), 2)
DTD = DTS // 86400
DTH = DTS // 3600 - DTD * 24
DTM = DTS // 60 - DTH * 60 - DTD * 24 * 60
DTS = np.round(DTS - (DTM + DTH * 60 + DTD * 24 * 60) * 60, 2)

print('\n')
print('\n')
print('**********************************************************************************************************')
print('**********************************************************************************************************')
print('**********************************************************************************************************')
print('\n')
print("%s days %s hours %s minutes and %s seconds overall required for full run" % (DTD, DTH, DTM, DTS))
print('\n')
print('\n')
