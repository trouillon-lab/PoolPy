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
import cmcrameri.cm as cmc




def plot_data(dir_WAs, max_diff, min_diff, start=0, stop=np.inf, step=1, x_axis='both', y_axis='all', **kwargs):
    N=start
    ls_names_met=['Method', 'Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps', 'N', 'diff' ]
    
    full_df_met=pd.DataFrame(columns=ls_names_met)

    while N<stop:
        Npath=os.path.join(dir_WAs,'N_'+str(N))
        if not os.path.exists(Npath):
            continue
        diff=min_diff
        while diff<=max_diff:
            start_time = time.time()
            dpath=os.path.join(Npath,'diff_'+str(diff))
            if not os.path.exists(dpath):
                continue
            metname=os.path.join(dpath, 'Metrics_N_'+str(N)+'_diff_'+str(diff)+'.csv')
            df_met=pd.read_csv(metname, header=True, index=True)
            df_met['N']=N
            df_met['diff']=diff

            diff+=1


        N+=step
                

