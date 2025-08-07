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
import json


from Functions import *
   



parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--step')
parser.add_argument('--dir_scramblers')
parser.add_argument('--dir_WAs')
parser.add_argument('--max_diff')
parser.add_argument('--max_dims')
parser.add_argument('--timeit')
parser.add_argument('--keep_ratios_constant')






args = parser.parse_args()



differentiate= 3 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
step= 10 if type(args.step)==type(None) else int(args.step)
keep_ratios_constant= False if type(args.keep_ratios_constant)==type(None) else args.keep_ratios_constant=='True' 
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'
max_diff= 4 if type(args.max_diff)==type(None) else int(args.max_diff)
max_dims= np.inf if type(args.max_dims)==type(None) else int(args.max_dims)



dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,
             'keep_ratios_constant': keep_ratios_constant,
             'start':start, 'stop':stop,  'step':step, 
             'dir_WAs':args.dir_WAs, 'dir_scramblers':args.dir_scramblers, 
             'max_differentiate': max_diff, 'max_dims':max_dims}


decode_sweep(**dict_kwargs)
