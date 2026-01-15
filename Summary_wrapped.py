import numpy as np
import argparse
import re
import itertools
import math
import pandas as pd
import os
import time

from Functions_wrapped import *


parser = argparse.ArgumentParser()
parser.add_argument('--start', type=int, default=10)
parser.add_argument('--stop', type=int, default=100)
parser.add_argument('--step', type=int, default=5)
parser.add_argument('--start_diff', type=int, default=1)
parser.add_argument('--stop_diff', type=int, default=10)
parser.add_argument('--step_diff', type=int, default=1)
parser.add_argument('--directory', type=str, default='./wrapped')
parser.add_argument('--folder_N', type=bool, default=True)
parser.add_argument('--folder_diff', type=bool, default=True)
parser.add_argument('--checks_hierarchical', type=int, default=65)
parser.add_argument('--max_prev', type=float, default=0.1)


args = parser.parse_args()
dict_kwargs = {
    'start': args.start,
    'stop': args.stop,
    'step': args.step,
    'start_diff': args.start_diff,
    'stop_diff': args.stop_diff,
    'step_diff': args.step_diff,
    'directory': args.directory,
    'folder_N': args.folder_N,
    'folder_diff': args.folder_diff,
    'checks_hierarchical': args.checks_hierarchical,
    'max_prev':args.max_prev
}

sweep_fly_summary(**dict_kwargs)
