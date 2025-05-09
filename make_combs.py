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


def iterative_add_N(dict_start, N_add, save=True,save_dir='./combinations/',
                     return_last=True, differentiate=3):
    tmp_d=copy.deepcopy(dict_start)
    N_start=dict_start[1][-1]
    i=0
    while i<N_add:
        print(N_start+i+1)
        tmp_d=add_1(tmp_d, ND=differentiate)
        if save:
            NM=os.path.join(save_dir,'N_'+str(N_start+i+1)+'.pk')
            with open(NM, 'wb') as handle:
                pickle.dump(tmp_d, handle, protocol=pickle.HIGHEST_PROTOCOL)
        i+=1
    if return_last:
        return(tmp_d)


parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--save_dir')


args = parser.parse_args()

f1n=os.path.join(args.save_dir,'N_'+str(args.start)+'.pk')

if os.path.isfile(f1n):
    with open(f1n, "rb") as input_file:
        dct_cmbn = pickle.load(input_file)

else:
    N=int(args.start)
    dct_cmbn={}
    dct_cmbn.update({1:np.arange(N)})
    for j in range(2,int(args.differentiate)+1):
        dct_cmbn.update({j:np.array(list(itertools.combinations(np.arange(N),j)))})

iterative_add_N(dict_start=dct_cmbn, N_add=int(args.stop)-int(args.start), 
                save_dir=args.save_dir, return_last=False, 
                differentiate=int(args.differentiate))









