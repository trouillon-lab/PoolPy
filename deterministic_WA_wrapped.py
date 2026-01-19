import numpy as np
import argparse
import time
import os

from Functions_wrapped import *


parser = argparse.ArgumentParser()
parser.add_argument('--start', type=int, default=50)
parser.add_argument('--stop', type=int, default=110)
parser.add_argument('--step', type=int, default=10)
parser.add_argument('--directory', type=str, default='./wrapped')
parser.add_argument('--max_diff', type=int, default=10)
parser.add_argument('--max_dims', type=int, default=np.inf)
parser.add_argument('--timeit', type=str, default='True')
parser.add_argument('--cleanup', type=str, default='True')
parser.add_argument('--max_prev', type=float, default=0.1)

args = parser.parse_args()

start_time = time.time()

start = args.start
stop = args.stop
step = args.step
save_dir = args.directory
max_diff = args.max_diff
max_dims = args.max_dims
timeit = args.timeit == 'True'
cleanup = args.cleanup == 'True'

dict_kwargs = {
    'differentiate': 2,
    'return_wa': True,
    'timeit': timeit,
    'start': start,
    'stop': stop,
    'step': step,
    'save_dir': save_dir,
    'max_diff': max_diff,
    'max_dims': max_dims,
    'max_prev':args.max_prev,
    'cleanup':cleanup
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
print("%s days %s hours %s minutes and %s seconds overall required for run" % (DTD, DTH, DTM, DTS))
print('\n')
print('\n')
