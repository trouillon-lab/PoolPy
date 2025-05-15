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


def is_consistent_precomp(well_assigner:np.array, differentiate:int, scrambler:dict) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    N=well_assigner.shape[0]
    for i in range(differentiate):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            this_sc=scrambler[diff]
            full_well_assigner=np.concatenate((full_well_assigner,np.bool(np.sum(well_assigner[this_sc], axis=1))))
    _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
    if len(counts)<full_well_assigner.shape[0]:
        return(False, full_well_assigner, counts)
    elif len(counts)==full_well_assigner.shape[0]:
        return(True,full_well_assigner, counts)
    else:
        print("Something is fishy")
        return(-1)
    
def is_consistent_precomp_alldiff(well_assigner:np.array, max_diff:int, scrambler:dict) -> list:
    if max_diff==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    N=well_assigner.shape[0]
    dict_counts={}
    for i in range(max_diff):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            this_sc=scrambler[diff]
            full_well_assigner=np.concatenate((full_well_assigner,np.bool(np.sum(well_assigner[this_sc], axis=1))))

        _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
        dict_counts.update({diff:counts})
    return dict_counts

def mean_metrics_precomp(well_assigner, differentiate, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent_precomp(well_assigner, differentiate, scrambler) 
    ET=extra_tests(counts)   
    rounds=np.sum(counts>1)/np.sum(counts>0)+1
    p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
    return BT+ET, ET,  rounds, p_check
    
#
def mean_metrics_precomp_alldiff(well_assigner, max_diff, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    dict_counts= is_consistent_precomp_alldiff(well_assigner, max_diff, scrambler) 
    dict_res={}
    for difo, counts in dict_counts.items():
        ET=extra_tests(counts)   
        rounds=np.sum(counts>1)/np.sum(counts>0)+1
        p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
        dict_res.update({difo:[BT+ET, ET,  rounds, p_check]})
    return dict_res

