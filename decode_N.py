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
parser.add_argument('--path_to_WA', type=str, default='./pooling_results', help="A string argument with default './pooling_results'")

args = parser.parse_args()

args_dict = vars(args)


