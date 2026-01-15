import numpy as np
import argparse
import re
import itertools
import math
import pandas as pd
import os
import time

def format_time(seconds):
    """Convert seconds to days, hours, minutes, seconds format."""
    days = int(seconds // 86400)
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    
    parts = []
    if days > 0:
        parts.append(f"{days}d")
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0:
        parts.append(f"{minutes}m")
    parts.append(f"{secs}s")
    
    return " ".join(parts)

def sweep_fly_summary(start, stop, step, start_diff, stop_diff, step_diff, directory, folder_N=False, folder_diff=False, checks_hierarchical=1000):
    """
    Calculate fly_summary for a range of N and differentiation values and save as CSV files.
    
    Parameters:
    -----------
    start : int
        Starting value for N
    stop : int
        Stopping value for N (inclusive)
    step : int
        Step size for N increments
    start_diff : int
        Starting value for differentiation
    stop_diff : int
        Stopping value for differentiation (inclusive)
    step_diff : int
        Step size for differentiation increments
    directory : str
        Base directory where files will be saved
    folder_N : bool, default=False
        If True, create a folder for each N value
    folder_diff : bool, default=False
        If True, create folders for each differentiation value
    checks_hierarchical : int, default=1000
        Number of checks for hierarchical method
    """
    
    # Create base directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)
    
    # Start overall timer
    start_time = time.time()
    
    # Loop over N values
    for N in range(start, stop + 1, step):
        N_start_time = time.time()
        
        # Determine N-level directory
        if folder_N:
            N_dir = os.path.join(directory, f'N_{N}')
            os.makedirs(N_dir, exist_ok=True)
        else:
            N_dir = directory
        
        # Loop over differentiation values
        for diff in range(start_diff, stop_diff + 1, step_diff):
            # Determine diff-level directory
            if folder_diff:
                if folder_N:
                    diff_dir = os.path.join(N_dir, f'diff_{diff}')
                else:
                    diff_dir = os.path.join(directory, f'diff_{diff}')
                os.makedirs(diff_dir, exist_ok=True)
            else:
                diff_dir = N_dir
            
            # Calculate metrics
            df = fly_summary(N, diff, checks_hierarchical=checks_hierarchical)
            
            # Save to CSV
            filename = f'Summary_N_{N}_diff_{diff}.csv'
            filepath = os.path.join(diff_dir, filename)
            df.to_csv(filepath, index=False)
        
        # Print time elapsed for this N
        N_elapsed = time.time() - N_start_time
        print(f"Completed N={N}, time elapsed {format_time(N_elapsed)}")
    
    # Print total completion time
    total_elapsed = time.time() - start_time
    print(f"\nAll calculations completed, total time elapsed {format_time(total_elapsed)}")

def IntegerToBinaryTF(num: int, ls_bn: list)-> list:
    if num >= 2:
        ls_bn=IntegerToBinaryTF(num // 2, ls_bn)
    ls_bn.append(num % 2==1)
    return(ls_bn)

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

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

def assign_wells_mat(n_compounds:int, **kwargs)->np.array:
    L1=np.ceil(np.sqrt(n_compounds))
    L2=L1-1 if L1*(L1-1)>=n_compounds else L1
    well_assigner=np.zeros((n_compounds, int(L1+L2)))==1
    for i in range(n_compounds):
        cp_id=[int(i//L2), int(L1+(i % L2))]
        well_assigner[i,cp_id]=True
    return(well_assigner)

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

def get_params_multidims(n_compounds:int, n_dims:int, **kwargs):
    L1=np.ceil(np.power(n_compounds, 1/n_dims))
    i=0
    while np.power(L1, n_dims-i-1)*np.power(L1-1, i+1)>=n_compounds:
        i=i+1
    ls_dim=[L1]*(n_dims-i)+[L1-1]*i
    return(ls_dim)

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
    WA=STD(N,q,k)
    return(WA)

def get_STD_params(n_compounds:int, differentiate=1, False_results=0, force_q=False, **kwargs):
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
    return(q,k)

def get_chinese_prameters(n_compounds:int,  differentiate:int, backtrack=False, special_diff=False, **kwargs)->np.array:
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
            npc=this_primes**this_exp
            if np.prod(npc)>=ND and np.sum(npc)<T:
                T=np.sum(npc)
                best_id=id_combo
        combo=ls_iter[best_id]
        carr=np.array(combo)
        flt=carr>0
        this_primes=nprimes[flt]
        this_exp=carr[flt]
        npc=this_primes**this_exp

        return npc
    
    if special_diff and differentiate==2:
        q=np.ceil(np.log(n_compounds)/np.log(3)).astype(int)
        t=int((q+5)*q/2)
        return q,t
    
    if special_diff and differentiate==3:
        q=np.ceil(np.log(n_compounds)/np.log(2)).astype(int)
        t=int((q-1)*q*2)
        return q,t

    return primes

def assign_wells_chinese_old(n_compounds:int,  differentiate:int,**kwargs)->np.array:
    prod=1
    n=1
    primes=[]
    c_id=np.arange(n_compounds) 
    while prod<n_compounds**differentiate:
        n=n+1
        if isprime(n):
            prod=prod*n
            primes.append(n)
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
            npc=this_primes**this_exp
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
                    WA[3*i+ii,j]=True if int(ls_nc3[j][i])==int(ii) else False
        k=3*q
        for i in range(q):
            for ii in range(i+1,q):
                for j in range(n_compounds):
                    WA[k,j]=True if int(ls_nc3[j][i])==int(ls_nc3[j][ii]) else False
                k+=1
        return q,t
    
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
                            WA[k,j]=True if int(ls_nc3[j][i])==nu and int(ls_nc3[j][ii])==nuu else False
                        k+=1
        return q,t


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

def generalized_factorial(N,pw):
    FT=1
    for i in range(N):
        FT*=(N-i)**(pw-1)
    return FT

def fly_summary(n_compounds:int,  differentiate:int, checks_hierarchical=65, **kwargs):
    """
    Calculate summary metrics for all pooling methods on the fly.
    Returns a DataFrame with one row per method.
    """
    methods = []

    
    # Multidimensional methods (2D, 3D, 4D)
    for n_dims in [2, 3, 4]:
        ls_dims=get_params_multidims(n_compounds, n_dims)
        arr_dims=np.array(ls_dims)
        method_name = f'multidim-{n_dims}'
        PCF=np.ones_like(arr_dims).astype(float)
        if differentiate>=np.max(arr_dims):
            PCT=100
        elif differentiate==1:
            PCT=0
        else:
            for i in range(differentiate):
                NG=(arr_dims-i-1)
                NT=(n_compounds-i-1)
                PC=NG/NT
                PCF*=PC
            PCT=100*(1-np.sum(PCF))


        MEE=np.min([differentiate**n_dims, generalized_factorial(differentiate,  n_dims), n_compounds ])-1
        method_dict = {
            'Pooling strategy': method_name,
            'Mean experiments': np.sum(arr_dims)+MEE,
            'Max samples per pool': np.prod(arr_dims)/np.min(arr_dims),
            'N pools': np.sum(arr_dims),
            'Percentage check': PCT,  
            'Mean extra experiments': MEE,
            'Mean steps': 1+PCT/100
        }
        methods.append(method_dict)
        
        # Matrix method - copy of 2D
        if n_dims == 2:
            matrix_dict = method_dict.copy()
            matrix_dict['Pooling strategy'] = 'Matrix'
            methods.append(matrix_dict)
    
    # Binary method
    method_name = 'Binary'
    if differentiate>=2:
        PC=100
        MEE=n_compounds
    else:
        PC=0
        MEE=0
    methods.append({
        'Pooling strategy': method_name,
        'Mean experiments': int(np.ceil(np.log2(n_compounds)))+MEE,  
        'Max samples per pool': int(np.ceil(n_compounds/2)),
        'N pools': int(np.ceil(np.log2(n_compounds))),
        'Percentage check': PC,  
        'Mean extra experiments': MEE,
        'Mean steps': 1+PC/100
    })
    
    # STD method
    method_name = 'STD'
    q,k = get_STD_params(n_compounds, differentiate)
    methods.append({
        'Pooling strategy': method_name,
        'Mean experiments': int(q*k),  
        'Max samples per pool': int(np.ceil(n_compounds/q)),
        'N pools':int(q*k),
        'Percentage check': 0,  
        'Mean extra experiments': 0,
        'Mean steps': 1
    })
    
    # Chinese Remainder method
    method_name = 'Chinese remainder'
    primez = get_chinese_prameters(n_compounds, differentiate)
    methods.append({
        'Pooling strategy': method_name,
        'Mean experiments': np.sum(primez),  
        'Max samples per pool': n_compounds/np.min(primez),
        'N pools': np.sum(primez),
        'Percentage check': 0,  
        'Mean extra experiments': 0,
        'Mean steps': 1
    })
    
    # Chinese Remainder backtrack method
    method_name = 'Ch. rm. bktrk'
    primez = get_chinese_prameters(n_compounds, differentiate, backtrack=True)
    methods.append({
        'Pooling strategy': method_name,
        'Mean experiments': np.sum(primez),  
        'Max samples per pool': n_compounds/np.min(primez),
        'N pools': np.sum(primez),
        'Percentage check': 0,  
        'Mean extra experiments': 0,
        'Mean steps': 1
    })
    
    # Chinese Remainder special method (only for differentiate 2 or 3)
    if differentiate in [2, 3]:
        method_name = 'Chinese special'
        q,t = get_chinese_prameters(n_compounds, differentiate, special_diff=True)
        if differentiate==2:
            MS=3**q
        elif differentiate==3:
            MS=2**(q-2)
        methods.append({
            'Pooling strategy': method_name,
            'Mean experiments': t,  
            'Max samples per pool':MS,
            'N pools': t,
            'Percentage check': 0,  
            'Mean extra experiments': 0 ,
            'Mean steps':1 
        })
        methods.append({
            'Pooling strategy': method_name,
            'Mean experiments': t,  
            'Max samples per pool':MS,
            'N pools': t,
            'Percentage check': 0,  
            'Mean extra experiments': 0,
            'Mean steps':1 
        })
        

    method_name = 'Hierarchical'
    metricas=calculate_metrics_hierarchical_fast(n_compounds, differentiate, checks=checks_hierarchical)
    methods.append({
            'Pooling strategy': method_name,
            'Mean experiments': metricas[0],  
            'Max samples per pool':metricas[1],
            'N pools': metricas[2],
            'Percentage check': metricas[3],  
            'Mean extra experiments': np.round(metricas[4],2),
            'Mean steps': metricas[5]
        })
    
    
    df = pd.DataFrame(methods)
    # Round all numeric columns to 2 decimal places
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df[numeric_cols] = df[numeric_cols].round(2)
    return df

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def pick_rand_pos(n_compounds,differentiate,max_checks=1e3):
    
    MC=0
    MP=1
    mp=0
    scaler=1
    ls_combs=[]
    ls_diffs=[]
    difo=differentiate
    while (MC<max_checks or MP>mp) and difo>0:
        MCI=math.comb(n_compounds,differentiate)
        MC+=MCI
        mp=MCI/MC
        ls_combs.append(MCI)
        ls_diffs.append(difo)
        difo-=1
    ls_posi=[]
    probi=np.array(ls_combs, dtype=float)
    probi/=np.sum(probi)
    differis=np.random.choice(ls_diffs, int(max_checks*scaler), p=probi)
    for differo in differis:
            rnd_pos=np.random.choice(np.arange(n_compounds), differo, replace=False)
            ls_posi.append(rnd_pos)

    return ls_posi


def iterative_splitter(id_samps, id_positives, ratio):
    
    if len(id_samps)<=ratio:
        return(len(id_samps))
    
    pools=list(split(id_samps, ratio))
    partials=0
    for pool in pools:
        if len(set(pool).intersection(id_positives))>0:
            partials+=iterative_splitter(pool,id_positives,ratio)
    return(ratio+partials)

def uneven_wrapper(n_samps, differentiate=-1):
    if differentiate==-1:
        differentiate=np.floor(n_samps/2)

    list_of_lists=[]
    def uneven_splits_maker(n_samps, previous_l):
        if n_samps<2:
            pass
        if len(previous_l)==0:
            C1=n_samps
        else:
            C1=previous_l[-1]
        MS=np.min([differentiate*5, C1 , n_samps])
        rationi=np.arange(2,MS+1)
        invariationi=n_samps/rationi
        linvariationi=np.array(list(set(invariationi.astype(int))))
        variationi=n_samps/linvariationi
        ratios=[rationi[np.argmin(np.abs(rationi-i))] for i in variationi ]
        for ratio in ratios:
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



def calculate_metrics_hierarchical_fast(n_compounds,  differentiate:int, checks=1e4, keep_ratios_constant=False,  **kwargs):
    id_samps=np.arange(n_compounds)
    details={}
    posiz=pick_rand_pos(n_compounds, differentiate, checks)
    if keep_ratios_constant:
        BM=[0,np.inf]
        for ratiof in np.arange(2,np.ceil(np.sqrt(n_compounds))):
            ratio=int(ratiof)
            NP=0
            FM=0
            for id_pos in posiz:
                posx=np.array(id_pos).reshape(-1,1)
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
            list_splits=uneven_wrapper(n_compounds, differentiate)
        ls_id=0
        for splito in list_splits:
            NP=0
            FM=0
            for id_pos in posiz:
                posx=np.array(id_pos)
                measures=iterative_uneven_splitter(id_samps,posx,splito)
                FM+=measures
                NP+=1
                    
            layers=len(splito)+1
            MC=int(np.ceil(n_compounds/splito[0]))
            #details.update({ls_id:[FM/NP, MC, splito, int(np.round((NP-1)/(NP),2)*100), FM/NP-splito[0],layers]})
            ls_id+=1
            if FM/NP<BM[1]:
                BM=[splito,FM/NP]
        layers=len(BM[0])+1
        MC=int(np.ceil(n_compounds/BM[0][0]))
        return([BM[1], MC, BM[0], int(np.round((NP-1)/(NP),2)*100), BM[1]-BM[0][0],layers])
    




parser = argparse.ArgumentParser()
parser.add_argument('--start', type=int, default=10)
parser.add_argument('--stop', type=int, default=100)
parser.add_argument('--step', type=int, default=5)
parser.add_argument('--start_diff', type=int, default=1)
parser.add_argument('--stop_diff', type=int, default=5)
parser.add_argument('--step_diff', type=int, default=1)
parser.add_argument('--directory', type=str, default='.')
parser.add_argument('--folder_N', type=bool, default=False)
parser.add_argument('--folder_diff', type=bool, default=False)
parser.add_argument('--checks_hierarchical', type=int, default=65)







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
    'checks_hierarchical': args.checks_hierarchical
}

sweep_fly_summary(**dict_kwargs)
