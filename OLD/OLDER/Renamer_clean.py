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
from Fast_functions import *
import argparse

parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--path', type=str, default='./', help='An integer argument with default 2')

args = parser.parse_args()
root_folder=args['path']


for dirpath, dirnames, filenames in os.walk(root_folder):
    # Exclude the safe folder if specified
    for filename in filenames:
        if filename.startswith('Random'):
