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
from Fast_functions import *



parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=2, help='An integer argument with default 2')
parser.add_argument('--n_samp', type=int, default=50, help='An integer argument with default 50')
parser.add_argument('--guesses', type=int, default=5, help='An integer argument for guesses of random WA with default 20')
parser.add_argument('--method', type=str, default='all', help="A string argument with default 'all'")
parser.add_argument('--path', type=str, default='./pooling_results', help="A string argument with default './pooling_results'")

args = parser.parse_args()

args_dict = vars(args)
