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
    while diff<(ND):
        new_part=np.vstack([combinantions_dictionary[diff].T,np.array([N]*len(combinantions_dictionary[diff]))])
        new_in=np.hstack([combinantions_dictionary[diff+1].T,new_part]).T
        new_cd.update({(diff+1):new_in})

        diff+=1
    

    return(new_cd)


def iterative_add_N(dict_start, N_add, save=True,save_dir='./combinations/',
                     return_last=True, differentiate=3, **kwargs):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    tmp_d=copy.deepcopy(dict_start)
    N_start=dict_start[1][-1]
    i=0
    while i<N_add:
        diri=os.path.join(save_dir,'N_'+str(N_start+i+2))
        if not os.path.exists(diri):
            os.makedirs(diri)
  
        print(N_start+i+2)
        tmp_d=add_1(tmp_d, ND=differentiate)
        if save:
            for ii in range(2,differentiate+1):
                this_diri=os.path.join(diri,'N_'+str(N_start+i+2)+'_diff_'+str(ii)+'.npz')
                this_diff=tmp_d[ii]
                #print(diri)
                np.savez_compressed(file=this_diri,ii=this_diff, allow_pickle=False)
        i+=1
    if return_last:
        return(tmp_d)


parser = argparse.ArgumentParser()
parser.add_argument('--differentiate')
parser.add_argument('--start')
parser.add_argument('--stop')
parser.add_argument('--save_dir')
parser.add_argument('--timeit')


args = parser.parse_args()


differentiate= 2 if type(args.differentiate)==type(None) else int(args.differentiate)
start= 50 if type(args.start)==type(None) else int(args.start)
stop= 110 if type(args.stop)==type(None) else int(args.stop)
save_dir= os.path.join(os.getcwd(),'outs') if type(args.save_dir)==type(None) else str(args.save_dir)
timeit= True if type(args.timeit)==type(None) else args.timeit=='True'


dict_kwargs={'differentiate':differentiate, 'return_wa':True, 'timeit':timeit,
             'start':start, 'stop':stop, 'save_dir':save_dir,'N_add':stop-start+1}


f1n=os.path.join(args.save_dir,'N_'+str(args.start)+'.pk')

if os.path.isfile(f1n) and False:
    with open(f1n, "rb") as input_file:
        dct_cmbn = pickle.load(input_file)

else:
    N=dict_kwargs['start']
    diff=dict_kwargs['differentiate']
    dct_cmbn={}
    dct_cmbn.update({1:np.arange(N)})
    for j in range(2,diff+1):
        dct_cmbn.update({j:np.array(list(itertools.combinations(np.arange(N),j)))})

iterative_add_N(dict_start=dct_cmbn, return_last=False, **dict_kwargs)









