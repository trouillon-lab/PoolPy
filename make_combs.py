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



def add_1(combinantions_dictionary, ND=5):
    N=combinantions_dictionary[1][-1]+1
    new_cd={1:np.append(combinantions_dictionary[1],N)}
    diff=1
    while diff<(ND-1):
        new_part=np.vstack([combinantions_dictionary[diff].T,np.array([N]*len(combinantions_dictionary[diff]))])
        new_in=np.hstack([combinantions_dictionary[diff+1].T,new_part]).T
        new_cd.update({(diff+1):new_in})

        diff+=1
    

    return(new_cd)


def iterative_add_N(dict_start, N_add, save=True,save_dir='./combinations/'):
    tmp_d=copy.deepcopy(dict_start)
    N_start=dict_start[1][-1]
    i=0
    while i<N_add:
        print(i)
        tmp_d=add_1(tmp_d)
        if save:
            NM=save_dir+'N_'+str(N_start+i+1)+'.pk'
            with open(NM, 'wb') as handle:
                pickle.dump(tmp_d, handle, protocol=pickle.HIGHEST_PROTOCOL)
        i+=1

    return(tmp_d)






