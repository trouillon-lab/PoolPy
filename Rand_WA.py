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


from Functions import *



def get_max_C(n_compounds, max_compounds):
    return int(n_compounds/2) if max_compounds==0 else max_compounds
def get_min_C(n_compounds, MC):
    return int(np.sqrt(n_compounds)) if int(np.sqrt(n_compounds))<MC else int(MC/2)

def get_max_W(n_compounds):
    return int(np.log2(n_compounds))
def get_min_W(n_compounds):
    return int(2*np.sqrt(n_compounds))