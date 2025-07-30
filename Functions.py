import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import pickle
import time
import os
import string

def int_to_base(n, N):
    """ Return base N representation for int n. """
    base_n_digits = string.digits + string.ascii_lowercase + string.ascii_uppercase
    result = ""
    if n < 0:
        sign = "-"
        n = -n
    else:
        sign = ""
    while n > 0:
        q, r = divmod(n, N)
        result += base_n_digits[r]
        n = q
    if result == "":
        result = "0"
    return sign + "".join(reversed(result))


#Coumpound counter starts from 1
# Helper function for binary translation
def IntegerToBinaryTF(num: int, ls_bn: list)-> list:
    if num >= 2:
        ls_bn=IntegerToBinaryTF(num // 2, ls_bn)
    ls_bn.append(num % 2==1)
    return(ls_bn)

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

# Function to assign each compound to its wells
def well_selecter(compound: int, n_wells:int, differentiate=1) -> np.array:
    if differentiate not in [1,2]:
        print('For the moment this code is only able to create well assignments matrices to distinguish up to combinations of 2 active compounds')
        return(-1)
    if differentiate==1:
        ls_bn=[]
        used_wells=IntegerToBinaryTF(compound, ls_bn)
        sel_wells=[False]*(n_wells-len(used_wells))+used_wells
    if differentiate==2:
        for i in range((n_wells-1)//3+1):
            if 0<compound and compound <= n_wells-1-3*i:
                sel_wells=(n_wells-1-3*i-compound)*[False]+[True]+i*3*[False]+[True]+[False]*(compound-1)
                break
            compound=compound-(n_wells-1-3*i)
    return(np.array(sel_wells))
    
def get_ncomp_from_nwells(nwells: int, differentiatie=1) ->int:
    if differentiatie==1:
        return(2**nwells-1)
    if differentiatie==2:
        temp_1=((nwells-2)//3)+1
        temp_2=temp_1-1
        return(nwells*temp_1-(temp_2*(temp_2+1)/2*3+temp_1)) #Some algebra gives this formula


# Functions to be called by user to create the compound-wells assignment matrix

''' Method 1: Pooling using matrix design'''
def assign_wells_mat(n_compounds:int, **kwargs)->np.array:
    L1=np.ceil(np.sqrt(n_compounds))
    L2=L1-1 if L1*(L1-1)>=n_compounds else L1
    well_assigner=np.zeros((n_compounds, int(L1+L2)))==1
    for i in range(n_compounds):
        cp_id=[int(i//L2), int(L1+(i % L2))]
        well_assigner[i,cp_id]=True
    return(well_assigner)


''' Method 2: Pooling using binary design'''
# This functions also identifies the minimum number of wells needed for the compounds and level of detail (differentiate) selected
def assign_wells_L(n_compounds:int, differentiate=1, **kwargs) -> np.array:

    if differentiate==1:
        n_wells=int(np.ceil(np.log2(n_compounds +1)))
    if differentiate==2:
        tentative=int(np.floor(np.sqrt(6*n_compounds))) #empirical evidence suggest almost this scaling, mathematical proof might arrive later
        for NW in [tentative-1,tentative,tentative+1]:
            if get_ncomp_from_nwells(NW, differentiatie=2)>=n_compounds:
                n_wells=NW
                break
        
    well_assigner=np.zeros((n_compounds, n_wells))==1
    for i in range(n_compounds):
        well_assigner[i,:]=well_selecter(i+1, n_wells, differentiate)
    return(well_assigner)

''' Method 2: Pooling using binary design'''
# This functions also identifies the minimum number of wells needed for the compounds and level of detail (differentiate) selected
def assign_wells_bin(n_compounds:int, differentiate=1, **kwargs) -> np.array:

    n_wells=int(np.ceil(np.log2(n_compounds +1)))
        
    well_assigner=np.zeros((n_compounds, n_wells))==1
    for i in range(n_compounds):
        well_assigner[i,:]=well_selecter(i+1, n_wells, differentiate=1)
    return(well_assigner)


''' Method 3: Pooling using multidimensional matrix design'''
def assign_wells_multidim(n_compounds:int, n_dims:int, **kwargs)->np.array:
    L1=np.ceil(np.power(n_compounds, 1/n_dims))
    i=0
    while np.power(L1, n_dims-i-1)*np.power(L1-1, i+1)>=n_compounds:
        i=i+1
    ls_dim=[L1]*(n_dims-i)+[L1-1]*i
    up_samps=np.prod(np.array(ls_dim))
    well_assigner=np.zeros((n_compounds, int(L1*(n_dims-i)+(L1-1)*i)))==1
    for j in range(n_compounds):
        cp_id=[]
        jj=np.copy(j)
        rem_dim=up_samps
        past_dims=0
        for k in range(n_dims):
            rem_dim=rem_dim/ls_dim[k]
            js=jj//rem_dim
            jj=jj-js*rem_dim
            jd=js+past_dims
            cp_id.append(int(jd))
            past_dims=past_dims+ls_dim[k]
        well_assigner[j,cp_id]=True
    return(well_assigner)




''' Method 4: Pooling using random design'''
def assign_wells_random(n_compounds:int,  differentiate:int, n_compounds_per_well=0, n_wells=0, guesses=0, Evaluate=False, return_me=False, **kwargs)->np.array:
    if guesses==0:
        guesses=n_compounds
    min_tests=np.inf

    if n_compounds_per_well==0 or n_wells==0:
        _,_, min_tests, WA_rand=find_rand_params(n_compounds=n_compounds, differentiate=differentiate, 
                                 n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, guesses=guesses)
        if return_me:
            return WA_rand,  min_tests
        
        return WA_rand
        


    if Evaluate:
        second_axis=np.tile(np.arange(n_wells),n_compounds_per_well).reshape(n_compounds_per_well,-1)
        for i in range(guesses):
            idt=np.random.randint(0,n_compounds,size=(n_compounds_per_well,n_wells) )
            well_assigner=np.zeros((n_compounds,n_wells))==1
            well_assigner[idt, second_axis]=True
            if guesses==1:
                if return_me:
                    mean_exp, _, _, p_check= mean_metrics(well_assigner, differentiate)
                    return well_assigner, mean_exp
                return well_assigner
            mean_exp, _, _, p_check= mean_metrics(well_assigner, differentiate)
            if p_check<1:
                if return_me:
                    return well_assigner,  mean_exp
                return well_assigner
            elif mean_exp<min_tests: 
                best_wa=well_assigner.copy()
                min_tests=mean_exp

        if return_me:
            return best_wa,  min_tests
        
        return best_wa

    _,_, min_tests, WA_rand=find_rand_params(n_compounds=n_compounds, differentiate=differentiate, 
                                 n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, guesses=guesses)
    if return_me:
        return WA_rand,  min_tests
    
    return WA_rand


''' Method 5: Pooling using STD design'''
# Hopeful implementation of the Shifted Transversal Design, what an acronym!

# From http://membres-timc.imag.fr/Nicolas.Thierry-Mieg/pdfs/ntm_IMA2012.pdf
# Corollary: If there are at most t positive variables and at most E false positive and E false negative observations: 
# STD(n;q;k) is a solution, when choosing q prime such that t⋅ Г(q,n)+2⋅ E  < q,  and  k = t⋅ Г+2⋅ E+1

def isprime(n:int):
    return re.compile(r'^1?$|^(11+)\1+$').match('1' * n) is None

def get_Gamma(q,N):
    return(int(np.ceil(np.log(N)/np.log(q))-1))

def get_s(N,j,q):
    vec=np.arange(N)
    out_vec=vec.copy()
    Gamma=get_Gamma(q,N)
    for ct in range(Gamma):
        c=ct+1
        out_vec=out_vec+j**c*(vec//q**c)
    return(out_vec)

def STD(N,q,k):
    L=np.zeros((k,q,N))==1
    for j in range(k):
        s=get_s(N,j,q)
        s=s%q
        idc=np.arange(N)
        L[j,s,idc]=True

    L=L.reshape(-1,N)
    return(L.T)

def assign_wells_STD(n_compounds:int, differentiate=1, False_results=0, force_q=False, **kwargs):
    N=n_compounds
    t=differentiate
    E=False_results
    poss_q=[x for x in range(n_compounds) if isprime(x)]
    for q in poss_q:
        if t*get_Gamma(q,N)+2*E<q:
            break
    if isprime(force_q):
        if t*get_Gamma(force_q,N)+2*E<force_q:
            q=force_q 
    Gamma=get_Gamma(q,N)
    k=t*Gamma+2*E+1
    if kwargs['return_wa']==False:
        return([q*k, np.sum((np.arange(n_compounds)%q)==0)])
    else:
        WA=STD(N,q,k)
        return(WA)
    
# Method from IMPROVED COMBINATORIAL GROUP TESTING ALGORITHMSFOR REAL-WORLD PROBLEM SIZES
# section The Chinese remainder sieve

def assign_wells_chinese(n_compounds:int,  differentiate:int, backtrack=False, special_diff=False, **kwargs)->np.array:
    prod=1
    n=1
    primes=[]
    c_id=np.arange(n_compounds) 
    while prod<n_compounds**differentiate:
        n=n+1
        if isprime(n):
            prod=prod*n
            primes.append(n)

    if backtrack:
        T=np.inf
        nprimes=np.array(primes)
        ND=n_compounds**differentiate
        ls_of_ls=[]
        LMP=np.log(primes[-1])
        for pi in primes:
            LE=np.floor(LMP/np.log(pi)).astype(int)
            ls_of_ls.append(list(range(LE+1)))
        ls_iter=list(itertools.product(*ls_of_ls))
        for id_combo, combo in enumerate(ls_iter):
            carr=np.array(combo)
            flt=carr>0
            this_primes=nprimes[flt]
            this_exp=carr[flt]
            npc=np.prod(this_primes**this_exp)
            if np.prod(npc)>=ND and np.sum(npc)<T:
                T=np.sum(npc)
                best_id=id_combo
        combo=ls_iter[best_id]
        carr=np.array(combo)
        flt=carr>0
        this_primes=nprimes[flt]
        this_exp=carr[flt]
        npc=this_primes**this_exp

        WA=np.zeros((np.sum(npc), n_compounds))==1
        past_primes=0
        for prime in npc:
            temp_wa=np.zeros((prime, n_compounds))==1
            for x in range(prime):
                ids=c_id%prime==x    
                temp_wa[x, ids]=True
            WA[past_primes:past_primes+prime,:]=temp_wa
            past_primes=past_primes+prime

        return(WA.T)   
    
    if special_diff and differentiate==2:
        q=np.ceil(np.log(n_compounds)/np.log(3)).astype(int)
        t=int((q+5)*q/2)
        WA=np.zeros((t, n_compounds))==1
        ls_nc3=[list(i)[::-1] for i in [int_to_base(j,3).zfill(q) for j in range(n_compounds)]]
        for i in range(q):
            for ii in range(3):
                for j in range(n_compounds):
                    WA[3*i+ii,j]=True if ls_nc3[j][i]==ii else False
        k=3*q
        for i in range(q):
            for ii in range(i+1,q):
                for j in range(n_compounds):
                    WA[k,j]=True if ls_nc3[j][i]==ls_nc3[j][ii] else False
                k+=1
        return(WA.T)
    
    if special_diff and differentiate==3:
        q=np.ceil(np.log(n_compounds)/np.log(2)).astype(int)
        t=int((q-1)*q*2)
        WA=np.zeros((t, n_compounds))==1
        ls_nc3=[list(i)[::-1] for i in [int_to_base(j,2).zfill(q) for j in range(n_compounds)]]
        k=0
        for i in range(q):
            for ii in range(i+1,q):
                for nu in [0,1]:
                    for nuu in [0,1]:
                        for j in range(n_compounds):
                            WA[k,j]=True if ls_nc3[j][i]==nu and ls_nc3[j][ii]==nuu else False
                        k+=1
        return(WA.T)




        



    WA=np.zeros((np.sum(primes), n_compounds))==1
    past_primes=0
    for prime in primes:
        temp_wa=np.zeros((prime, n_compounds))==1
        for x in range(prime):
            ids=c_id%prime==x    
            temp_wa[x, ids]=True
        WA[past_primes:past_primes+prime,:]=temp_wa
        past_primes=past_primes+prime

    return(WA.T)

def iterative_splitter(id_samps, id_positives, ratio):
    
    if len(id_samps)<=ratio:
        return(len(id_samps))
    
    pools=list(split(id_samps, ratio))
    partials=0
    for pool in pools:
        if len(set(pool).intersection(id_positives))>0:
            partials+=iterative_splitter(pool,id_positives,ratio)
    return(ratio+partials)

def uneven_wrapper(n_samps):
    list_of_lists=[]
    def uneven_splits_maker(n_samps, previous_l):
        if n_samps<2:
            pass
        for ratio in np.arange(2,np.floor(n_samps/2)+1):
            this_l=previous_l.copy()
            this_l.append(int(ratio))
            list_of_lists.append(this_l)
            uneven_splits_maker(np.ceil(n_samps/ratio), this_l)

    uneven_splits_maker(n_samps,[])
    
    
    return(list_of_lists)
        

def iterative_uneven_splitter(id_samps, id_positives, ratios):
    if len(ratios)==1:
        ratio=ratios[0]
        ratios=[np.inf]
    else:
        ratio=ratios[0]
        ratios=ratios[1:]

    if len(id_samps)<=ratio:
        return(len(id_samps))

    pools=list(split(id_samps, ratio))
    partials=0
    for pool in pools:
        if len(set(pool).intersection(id_positives))>0:
            partials+=iterative_uneven_splitter(pool,id_positives,ratios)
    return(ratio+partials)







def calculate_metrics_hierarchical(n_compounds,  differentiate:int,  **kwargs):
    keep_ratios_constant=kwargs['keep_ratios_constant']
    id_samps=np.arange(n_compounds)
    details={}
    
    if keep_ratios_constant:
        BM=[0,np.inf]
        for ratiof in np.arange(2,np.ceil(np.sqrt(n_compounds))):
            ratio=int(ratiof)
            NP=0
            FM=0
            for n_pos in np.arange(differentiate+1):
                for id_pos in itertools.combinations(np.arange(n_compounds),n_pos):
                    posx=np.array(id_pos)
                    measures=iterative_splitter(id_samps,posx,ratio)
                    FM+=measures
                    NP+=1
                    
            layers=int(np.ceil(np.log(n_compounds)/np.log(ratio)))
            MC=int(np.ceil(n_compounds/ratio))

            details.update({ratio:[BM[1], MC, BM[0], int(np.round((NP-1)/(NP),2)*100), BM[1]-BM[0],layers]})
            if FM/NP<BM[1]:
                BM=[ratio,FM/NP]
        layers=int(np.ceil(np.log(n_compounds)/np.log(BM[0])))
        MC=int(np.ceil(n_compounds/BM[0]))
        return([BM[1], MC, BM[0], int(np.round((NP-1)/(NP),2)*100), BM[1]-BM[0],layers, details ]) 

    else:
        BM=[[0],np.inf]
        
        if 'ls_splits' in kwargs.keys():
            list_splits=[kwargs['ls_splits']]
        else:
            list_splits=uneven_wrapper(n_compounds)
        ls_id=0
        for splito in list_splits:
            NP=0
            FM=0
            for n_pos in np.arange(differentiate+1):
                for id_pos in itertools.combinations(np.arange(n_compounds),n_pos):
                    posx=np.array(id_pos)
                    measures=iterative_uneven_splitter(id_samps,posx,splito)
                    FM+=measures
                    NP+=1
                    
            layers=len(splito)+1
            MC=int(np.ceil(n_compounds/splito[0]))
            details.update({ls_id:[BM[1], MC, BM[0], int(np.round((NP-1)/(NP),2)*100), BM[1]-BM[0][0],layers]})
            ls_id+=1
            if FM/NP<BM[1]:
                BM=[splito,FM/NP]
        layers=len(BM[0])+1
        MC=int(np.ceil(n_compounds/BM[0][0]))
        return([BM[1], MC, BM[0], int(np.round((NP-1)/(NP),2)*100), BM[1]-BM[0][0],layers, details ])




        

''' Additional Functions '''

# Helper function to actually implement the experiment
def from_well_get_compuonds(well:int, well_assigner: np.array)->np.array :
    return(np.array(np.where(well_assigner[:,well-1]))[0]+1)

def from_compound_get_wells(compound: int, well_assigner: np.array)-> np.array:
    return(np.array(np.where(well_assigner[compound-1,:]))[1]+1)

# Consistency check

# def is_consistent_old(well_assigner:np.array, differentiate:int) -> list:
#     n_comp=well_assigner.shape[0]
#     if differentiate==1:
#         full_well_assigner=well_assigner.copy()
#     if differentiate==2:
#         n_perms=int(n_comp*(n_comp-1)/2+n_comp)
#         full_well_assigner=np.zeros((n_perms, well_assigner.shape[1]))==1
#         full_well_assigner[-n_comp:,:]=well_assigner.copy()
#         for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),2)):
#             k,l=j
#             full_well_assigner[i,:]=well_assigner[k,:]+well_assigner[l,:]
#     if differentiate==3:
#         n_perms=int(n_comp*(n_comp-1)*(n_comp-2)/6+n_comp*(n_comp-1)/2+n_comp)
#         n3=int(n_comp*(n_comp-1)*(n_comp-2)/6)
#         full_well_assigner=np.zeros((n_perms, well_assigner.shape[1]))==1
#         full_well_assigner[-n_comp:,:]=well_assigner.copy()
#         for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),3)):
#             k,l,m=j
#             full_well_assigner[i,:]=well_assigner[k,:]+well_assigner[l,:]+well_assigner[m,:]
#         for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),2)):
#             k,l=j
#             full_well_assigner[i+n3,:]=well_assigner[k,:]+well_assigner[l,:]
#     _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
#     if np.unique(full_well_assigner, axis=0).shape[0]<full_well_assigner.shape[0]:
#         return(False, full_well_assigner, counts)
#     elif np.unique(full_well_assigner, axis=0).shape[0]==full_well_assigner.shape[0]:
#         return(True,full_well_assigner, counts)
#     else:
#         return("Something is fishy")
    

#
def is_consistent(well_assigner:np.array, differentiate:int) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    N=well_assigner.shape[0]
    for i in range(differentiate):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            N_cmbn=math.comb(N,diff)
            temp_well_assigner=np.zeros((N_cmbn, well_assigner.shape[1]))==1
            for l,k in enumerate(itertools.combinations(np.arange(N),diff)):
                temp_well_assigner[l,:]=np.sum(well_assigner[k,:], axis=0).reshape(1,-1)
            full_well_assigner=np.concatenate((full_well_assigner,temp_well_assigner))
    _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
    if np.unique(full_well_assigner, axis=0).shape[0]<full_well_assigner.shape[0]:
        return(False, full_well_assigner, counts)
    elif np.unique(full_well_assigner, axis=0).shape[0]==full_well_assigner.shape[0]:
        return(True,full_well_assigner, counts)
    else:
        print("Something is fishy")
        return(-1)
    
#
def decode(well_assigner:np.ndarray, readout:np.ndarray, differentiate:int) -> list:
    N=well_assigner.shape[0]
    for i in range(differentiate):
        resulti=[]
        diff=i+1
        if diff ==1:
            resulti.extend(list(range(well_assigner.shape[0])))
            full_well_assigner=well_assigner.copy()
        else:
            N_cmbn=math.comb(N,diff)
            temp_well_assigner=np.zeros((N_cmbn, well_assigner.shape[1]))==1
            for l,k in enumerate(itertools.combinations(np.arange(N),diff)):
                resulti.append(k)
                temp_well_assigner[l,:]=np.sum(well_assigner[k,:], axis=0).reshape(1,-1)
            full_well_assigner=np.concatenate((full_well_assigner,temp_well_assigner))
    idxs=[i for i in range(full_well_assigner.shape[0]) if np.array_equal(full_well_assigner[i,:],readout)]
    if len(idxs)==0:
        print('No match')
        return(-1)
    
    return [resulti[i] for i in idxs]
        
def extra_tests(counts:np.array)->float:
    return(np.sum(counts*(counts-1))/np.sum(counts))

def extra_test_corrected(counts:np.array, N:int)->float:
    MC=np.array([N]*len(counts))
    max_c=np.minimum(counts,MC)
    return(np.sum(counts*(max_c)*((counts-1)>0))/np.sum(counts))
    
def mean_tests(well_assigner, differentiate, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent(well_assigner, differentiate)    
    ET=extra_tests(counts)
    return(BT+ET)


def mean_metrics(well_assigner, differentiate, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent(well_assigner, differentiate) 
    ET=extra_tests(counts)   
    rounds=np.sum(counts>1)/np.sum(counts>0)+1
    p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
    return BT+ET, ET,  rounds, p_check


def find_dims(n_compounds,differentiate, **kwargs):
    ND=2
    ndmin=2
    min_tests=n_compounds**2
    while 2**ND<n_compounds:
        WA_tmp=assign_wells_multidim(n_compounds,ND)
        mean_exp, _,  _, _= mean_metrics(WA_tmp, differentiate)
        if mean_exp<min_tests:
            ndmin=ND
            min_tests=mean_exp
        ND=ND+1
    return ndmin





def find_rand_params(n_compounds:int, differentiate:int, n_compounds_per_well=0, n_wells=0, guesses=0, 
                     max_compounds=0, max_redundancy=4, min_redundancy=1):
    skip_compounds=True
    skip_wells=True
    if n_compounds_per_well==0:
        skip_compounds=False
    if n_wells==0:
        skip_wells=False
    if guesses==0:
        guesses=n_compounds

    MC= int(n_compounds/2) if max_compounds==0 else max_compounds
    mc=int(np.sqrt(n_compounds)) if int(np.sqrt(n_compounds))<MC else int(MC/2)
    arr_comp=np.arange(int(mc),int(MC+1))
    mw=int(0.75*np.sqrt(n_compounds))
    MW=int(2*np.sqrt(n_compounds))
    while MW-mw<10:
        mw=int(abs(mw-1))
        MW=int(MW+1)

    arr_wells=np.arange(mw,MW)
    min_tests=np.inf
    for comp in arr_comp:
        if skip_compounds:
            comp=n_compounds_per_well

        for wells in arr_wells:
            if skip_wells:
                if skip_compounds:
                    
                    return n_compounds_per_well, n_wells, assign_wells_random(n_compounds=n_compounds, differentiate=differentiate, 
                                               n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, guesses=guesses, Evaluate=True)
                wells=n_wells
                
            if comp*wells>max_redundancy*n_compounds or comp*wells<min_redundancy*n_compounds: continue 
            WA_tmp, mean_exp=assign_wells_random(n_compounds, differentiate, comp, wells, guesses, Evaluate=True, return_me=True)
            if mean_exp<min_tests:
                Comp=comp
                Wells=wells
                min_tests=mean_exp
                min_wa=WA_tmp
            if skip_wells:
                break
        if skip_compounds:
            break

    return Comp, Wells, min_tests, min_wa
    
                       



# Function to export design table as .csv file
def output_table(well_assigner:np.array,output_file_name='output_table'):
    output_df=pd.DataFrame(well_assigner)
    '''
    compound_list=[f'compound_{x+1}' for x in output_df.index]
    pool_list=[f'pool_{x+1}' for x in output_df.columns]
    output_df.index=compound_list
    output_df.columns=pool_list
    '''
    output_df.rename(columns=lambda x: 'compound_'+str(x+1), index=lambda x: 'pool_'+str(x+1), inplace=True)

    output_df.to_csv(output_file_name+'.csv')
    print('design table saved as ',output_file_name,'.csv')
    
    return(output_df)


# Function to extracts metric values from a given file name generated with inline_print
def extract_metrics(filename):
    
    pattern = r"diff_(\d+)_NS_(\d+)_NW_(\d+)_MS_(\d+)"
    match = re.match(pattern, filename)
    if match:
        return {
            "diff": int(match.group(1)),
            "NS": int(match.group(2)),
            "NW": int(match.group(3)),
            "MS": int(match.group(4))
        }
    else:
        raise ValueError(f"Filename '{filename}' does not match the expected format.")



''' Method comparison '''

def full_method_comparison(**kwargs):
    methods=['matrix', 'random', 'STD', 'Chinese trick']
    # matrix assignment
    

    # multidimensional matrix
    if 'n_dims' in kwargs.keys():
        WA_mul=assign_wells_multidim(**kwargs)
        ndmin=kwargs['n_dims']
        multi=['multidim: '+str(ndmin)]
        WA_list=[WA_mul]
    elif 'all_dims' in kwargs.keys():
        if kwargs['all_dims']:
            WA_list=[]
            multi=[]
            for i in np.arange(2,int(np.ceil(np.log(kwargs['n_compounds'])/np.log(2)))):
                if i>kwargs['max_dims']:
                    continue
                WA_mul=assign_wells_multidim(n_dims=i, **kwargs)
                WA_list.append(WA_mul)
                multi.append('multidim: '+str(i))

        else:
            ndmin= find_dims(**kwargs)
            WA_mul=assign_wells_multidim(n_dims=ndmin, **kwargs)
            multi=['multidim: '+str(ndmin)]
            WA_list=[WA_mul]
            


    else:
        ndmin= find_dims(**kwargs)
        WA_mul=assign_wells_multidim(n_dims=ndmin, **kwargs)
        multi=['multidim: '+str(ndmin)]
        WA_list=[WA_mul]
    
    multi.extend(methods)
    methods=multi.copy()

    WA_mat=assign_wells_mat(**kwargs)

    # random assignment

    WA_ran=assign_wells_random(**kwargs)

    # STD asignment 
    WA_std=assign_wells_STD(**kwargs)

    

    # chinese trick assignment
    WA_chin=assign_wells_chinese(**kwargs)


    WA_list.extend([WA_mat, WA_ran,WA_std, WA_chin])

    if kwargs['differentiate']<2:
        WA_bin=assign_wells_L(**kwargs)
        methods.append('Binary')
        WA_list.append(WA_bin)

    #hierarchical
        
    Hier=calculate_metrics_hierarchical(**kwargs)
    #return([BM[0], BM[1],layers, MC, details ])

    ls_met=[]
    ls_names_met=['mean_experiments', 'max_compounds_per_well', 'n_wells', 'percentage_check', 'mean_extra_exp', 'mean_steps']
    for method, WA in zip(methods, WA_list):
        mean_exp, extra_exp,  _, perc_check= mean_metrics(WA, **kwargs)
        n_wells=WA.shape[1]
        M_exp=np.round(mean_exp, 2)
        max_comp=np.max(np.sum(WA, axis=0))
        ls_met.append([M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100])
    ls_met.append(Hier[:-1])
    full_methods=methods.copy()
    full_methods.append('Hierarchical')
    df_met=pd.DataFrame(ls_met)

    dict_wa={method: WA for method, WA in zip(methods, WA_list)}
    dict_wa.update({'Hierarchical':Hier[5]})

    idx_renamer={i:j for i,j in zip(df_met.index, full_methods)}
    col_renamer={i:j for i,j in zip(df_met.columns, ls_names_met)}
    df_met.rename(index=idx_renamer, columns=col_renamer, inplace=True)

    ret_wa= kwargs['return_wa'] 
    if ret_wa:
        return df_met, dict_wa
    return df_met

def full_sweep_comparison(start=50, stop=150, step=10, **kwargs):
    dict_comp={}
    current=start
    kwargs['return_wa']=True
    while current<stop:
        time0=time.time()
        if kwargs['timeit']:
            print(current)
        df_met, dict_wa=full_method_comparison(n_compounds=current, **kwargs)
        dict_comp.update({current:[df_met, dict_wa]})
        current=current+step
        if kwargs['timeit']:
            print("segment time: %s seconds" % np.round(time.time() - time0, 1))
    return dict_comp

def single_method_sweep(start=50, stop=150, step=10, **kwargs):
    dict_comp={'metrics_reading_key':['mean_experiments', 'max_compounds_per_well', 'n_wells', 
                                      'percentage_check', 'mean_extra_exp', 'mean_steps']}
    if kwargs['inline_print']:
        fpath=os.path.join(kwargs['base_dir'],kwargs['method'])
        if not os.path.exists(fpath):
            os.makedirs(fpath)
    current=start
    while current<stop:
        time0=time.time()
        if kwargs['timeit']:
            print(current)
        match kwargs['method']:
            case 'STD':
                WA=assign_wells_STD(n_compounds=current, **kwargs)
                if kwargs['return_wa']:
                    mean_exp, extra_exp,  _, perc_check= mean_metrics(WA, **kwargs)
                    n_wells=WA.shape[1]
                    M_exp=np.round(mean_exp, 2)
                    max_comp=np.max(np.sum(WA, axis=0))
                    dict_wa={'WA': WA, 'metrics':[M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100]}

                elif kwargs['inline_print']:
                    full_file_dir=os.path.join(fpath,'diff_'+str(kwargs['differentiate'])+'_NS_'+
                                               str(current)+'_NW_'+str(WA[0])+'_MS_'+str(WA[1])+".txt")
                    open(full_file_dir, 'a').close()
                else:
                    dict_wa={'metrics':[WA[0], WA[1], WA[0], 0,  0,]}
            case 'CT':
                WA=assign_wells_chinese(n_compounds=current, **kwargs)
                if kwargs['return_wa']:
                    mean_exp, extra_exp,  _, perc_check= mean_metrics(WA,  **kwargs)
                    n_wells=WA.shape[1]
                    M_exp=np.round(mean_exp, 2)
                    max_comp=np.max(np.sum(WA, axis=0))
                    dict_wa={'WA': WA, 'metrics':[M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100]}
                elif kwargs['inline_print']:
                    full_file_dir=os.path.join(fpath,'diff_'+str(kwargs['differentiate'])+'_NS_'+str(current)+
                                               '_NW_'+str(WA)+ '_MS_'+str(int(np.ceil(current/2)))+".txt")
                    open(full_file_dir, 'a').close()
                else:
                    dict_wa={'metrics':[WA, int(np.ceil(current/2)), WA, 0,  0,]}

            case 'random':
                WA=assign_wells_random(n_compounds=current, **kwargs)
                mean_exp, extra_exp,  _, perc_check= mean_metrics(WA, **kwargs)
                n_wells=WA.shape[1]
                M_exp=np.round(mean_exp, 2)
                max_comp=np.max(np.sum(WA, axis=0))
                dict_wa={'WA': WA, 'metrics':[M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100]}

            case 'matrix':
                WA=assign_wells_mat(n_compounds=current, **kwargs)
                mean_exp, extra_exp,  _, perc_check= mean_metrics(WA, **kwargs)
                n_wells=WA.shape[1]
                M_exp=np.round(mean_exp, 2)
                max_comp=np.max(np.sum(WA, axis=0))
                dict_wa={'WA': WA, 'metrics':[M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100]}

            case 'multidim':
                if 'n_dims' in kwargs.keys():
                    WA=assign_wells_multidim(n_compounds=current, **kwargs)
                    ndmin=kwargs['n_dims']
                else:
                    ndmin= find_dims(**kwargs)
                    WA=assign_wells_multidim(n_dims=ndmin, n_compounds=current,**kwargs)
                mean_exp, extra_exp,  _, perc_check= mean_metrics(WA, **kwargs)
                n_wells=WA.shape[1]
                M_exp=np.round(mean_exp, 2)
                max_comp=np.max(np.sum(WA, axis=0))
                dict_wa={'WA': WA, 'metrics':[M_exp, max_comp, n_wells, int(perc_check),  extra_exp,1+perc_check/100]}

            case 'hierarchical':
                Hier=calculate_metrics_hierarchical(n_compounds=current, **kwargs)
                dict_wa={'WA': Hier[5], 'metrics':Hier[:-1]}



        
        if kwargs['timeit']:
            print("segment time: %s seconds" % np.round(time.time() - time0, 1))
        if not kwargs['inline_print'] or kwargs['method'] not in ['STD', 'CT']:
            dict_comp.update({current:dict_wa})

        current=current+step


    return(dict_comp)

    




######################################################################################################
######################################################################################################
###########################                                                ###########################
###########################       BEGINNING OF PREVALENCE SECTION          ###########################
###########################                                                ###########################
######################################################################################################
######################################################################################################







def get_diff_from_prev(n_compounds, prev, p_cut=1e-4):
    ls_p=[]
    for i in range(n_compounds):
        combinations=math.comb(n_compounds,i)
        P=prev**i*(1-prev)**(n_compounds-i)*combinations
        ls_p.append(P)
        if np.sum(ls_p) > 1- p_cut:
            break
        
    return(np.array(ls_p)/np.sum(ls_p))

def prevalence_pooling(n_compounds, prev,  p_cut=1e-4, file='Final_precomputed_file.pk', **kwargs):
    with open(file, 'rb') as handle:
        f1=pickle.load(handle)
    Ps=get_diff_from_prev(n_compounds, prev,  p_cut)
    max_diff=len(Ps)-1
    if max_diff > 4:
        sys.exit(f'Maximum number of positives ({max_diff}) needed for requested prevalnce too high. The precomputed maximum is 4.')

    else:
        with open(file, 'rb') as handle:
            f1=pickle.load(handle)
        f2=f1['Differentiate '+str(max_diff)]
        a2=np.array(list(f2))
        md=np.max(a2)
        if n_compounds>md:
            sys.exit(f'Maximum number of samples ({n_compounds}) too high for the chosen number of positives \n Below are displayed the information for {md} samples and up to {max_diff} positives. \n You can run the code following the commands below for your specific case.')
        elif np.sum(a2==n_compounds)==0:
            print(f'There is no precomputed strategy for {n_compounds} samples. \n The closest precomputed strategy is for {md} samples with up to {max_diff} positives')
            md=np.min(a2[a2>=n_compounds])
        else:
            md=n_compounds

    CR=f2[md]
    FCR=recalculate_full_metrics(CR, max_diff, Ps, n_compounds=n_compounds, **kwargs)
    return(FCR)

def recalculate_full_metrics(full_ls, differentiate, Ps, **kwargs):
    summary=full_ls[0]
    dict_WA=full_ls[1]
    ls_met=[]
    ls_names_met=['mean_experiments', 'max_compounds_per_well', 'percentage_check', 'mean_extra_exp', 'mean_steps']
    methods=[]
    dict_ps={}
    final_df=pd.DataFrame()
    if 'Hierarchical' in dict_WA.keys():
        ls_split=summary.loc['Hierarchical', 'n_wells']
    for diff in range (differentiate+1):
        methods=[]
        ls_met=[]
        NW=[]
        for i,j in dict_WA.items():
            WA=j
            if i=='Hierarchical':
                Hier=calculate_metrics_hierarchical(differentiate= diff, ls_split=[ls_split],**kwargs)
                H1=Hier[:2]
                H2=Hier[3:-1]
                H1.extend(H2)
                ls_met.append(H1)
                NW.append(Hier[2])
            else:
                mean_exp, extra_exp,  _, perc_check= mean_metrics(WA,differentiate= diff, **kwargs)
                n_wells=WA.shape[1]
                M_exp=np.round(mean_exp, 2)
                max_comp=np.max(np.sum(WA, axis=0))
                ls_met.append([M_exp, max_comp, perc_check,  extra_exp,1+perc_check/100])
                NW.append(n_wells)
            methods.append(i)

        df_met=pd.DataFrame(ls_met)
        idx_renamer={i:j for i,j in zip(df_met.index, methods)}
        col_renamer={i:j for i,j in zip(df_met.columns, ls_names_met)}
        df_met.rename(index=idx_renamer, columns=col_renamer, inplace=True)
        if final_df.shape==df_met.shape:
            final_df=final_df+(df_met*Ps[diff])
        else:
            final_df=(df_met*Ps[diff])

        dict_ps.update({diff:[Ps[diff],df_met]})
    final_df['n_wells']=NW
    dict_ps.update({diff:[Ps[diff],df_met]})
    dict_out={'Prevalence_scale':final_df}
    dict_out.update(dict_ps)
    dict_out.update({diff:[Ps[diff],df_met]})
    return(dict_out)
    
