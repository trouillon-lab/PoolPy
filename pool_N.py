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
import argparse

parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=2, help='An integer argument with default 2')
parser.add_argument('--N', type=int, default=50, help='An integer argument with default 50')
parser.add_argument('--method', type=str, default='all', help="A string argument with default 'all'")
parser.add_argument('--path', type=str, default='./', help="A string argument with default './'")

args = parser.parse_args()

args_dict = vars(args)

diff= args_dict['differentiate']
N=args_dict['N']


scrambler={1:np.arange(N)}
for j in range(2,diff+1):
    scrambler.update({j:np.array(list(itertools.combinations(np.arange(N),j)))})