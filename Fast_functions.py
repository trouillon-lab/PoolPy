import numpy as np
import re
import itertools
import math
import pandas as pd
import os
import time
import string
import shutil
import scipy
from sklearn.linear_model import ElasticNet
from scipy.optimize import minimize
from scipy.optimize import nnls
from sklearn.linear_model import Lasso

def int_to_base(n, N):
    """Return base N representation for int n."""
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

def str_to_tuple(string, delimiter='-'):
    return tuple(string.split(delimiter))

def tuple_to_str(tuple_type, delimiter='-'):
    return delimiter.join(map(str,tuple_type))

# Utility: Generate a sparse random array of length L with p nonzero elements (values in [0, 3))
def sparse_uniform_array(L, p, low=0.0, high=3.0, random_state=None):
    """
    Generate an array of length L with exactly p nonzero elements, each drawn uniformly from [low, high).
    The positions of the nonzero elements are chosen at random.
    """
    if p > L:
        raise ValueError("p cannot be greater than L")
    rng = np.random.default_rng(random_state)
    arr = np.zeros(L)
    if p > 0:
        idx = rng.choice(L, size=p, replace=False)
        arr[idx] = rng.uniform(low, high, size=p)
    return arr

def sweep_fly_summary(start, stop, step, start_diff, stop_diff, step_diff, directory, folder_N=False, folder_diff=False, checks_hierarchical=1000, max_prev=0.1, **kwargs):
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

    random_filename_re = re.compile(r"^Random_diff_(?P<diff>[0-9]+)_NS_(?P<ns>[0-9]+)_NW_(?P<nw>[0-9.]+)_MS_(?P<ms>[0-9.]+)_PC_(?P<pc>[0-9.]+)_ME_(?P<me>[0-9.]+)\.txt$")

    def parse_random_metrics(candidates, N, diff):
        for cand in candidates:
            if not os.path.isdir(cand):
                continue
            for fname in os.listdir(cand):
                match = random_filename_re.match(fname)
                if not match:
                    continue
                parts = match.groupdict()
                if int(parts['diff']) != diff or int(parts['ns']) != N:
                    continue
                nw = float(parts['nw'])
                ms = float(parts['ms'])
                pc = float(parts['pc'])
                me = float(parts['me'])
                extra = max(me - nw, 0.0)
                return {
                    'Pooling strategy': 'Random',
                    'Max experiments': me,
                    'Max samples per pool': ms,
                    'N pools': nw,
                    'Percentage check': pc,
                    'Max extra experiments': extra,
                    'Mean steps': 1 + pc/100.0
                }
        return None
    
    # Loop over N values
    for N in range(start, stop + 1, step):
        N_start_time = time.time()
        
        # Determine N-level directory
        if folder_N:
            N_dir = os.path.join(directory, f'N_{N}')
            os.makedirs(N_dir, exist_ok=True)
        else:
            N_dir = directory
        stop_diffo=stop_diff if stop_diff>0 else int(max_prev*N)+2
        # Loop over differentiation values
        for diff in range(start_diff, stop_diffo + 1, step_diff):
            if diff>max_prev*N:
                break
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

            # If a random design summary file exists, append its metrics
            candidate_dirs = list(dict.fromkeys([
                diff_dir,
                os.path.join(N_dir, f'diff_{diff}'),
                os.path.join(directory, f'N_{N}', f'diff_{diff}')
            ]))
            random_row = parse_random_metrics(candidate_dirs, N, diff)
            if random_row is not None:
                df = pd.concat([df, pd.DataFrame([random_row])], ignore_index=True)
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                df[numeric_cols] = df[numeric_cols].round(2)
            
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
        T=np.prod(primes)
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
        try:
            combo=ls_iter[best_id]
            carr=np.array(combo)
            flt=carr>0
            this_primes=nprimes[flt]
            this_exp=carr[flt]
            npc=this_primes**this_exp
        except:
            npc=nprimes


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
        try:
            combo=ls_iter[best_id]
            carr=np.array(combo)
            flt=carr>0
            this_primes=nprimes[flt]
            this_exp=carr[flt]
            npc=this_primes**this_exp
        except:
            npc=nprimes

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
        return (WA.T)
    
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
        return (WA.T)


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
            'Max experiments': np.sum(arr_dims)+MEE,
            'Max samples per pool': np.prod(arr_dims)/np.min(arr_dims),
            'N pools': np.sum(arr_dims),
            'Percentage check': PCT,  
            'Max extra experiments': MEE,
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
        'Max experiments': int(np.ceil(np.log2(n_compounds)))+MEE,  
        'Max samples per pool': int(np.ceil(n_compounds/2)),
        'N pools': int(np.ceil(np.log2(n_compounds))),
        'Percentage check': PC,  
        'Max extra experiments': MEE,
        'Mean steps': 1+PC/100
    })
    
    # STD method
    method_name = 'STD'
    q,k = get_STD_params(n_compounds, differentiate)
    methods.append({
        'Pooling strategy': method_name,
        'Max experiments': int(q*k),  
        'Max samples per pool': int(np.ceil(n_compounds/q)),
        'N pools':int(q*k),
        'Percentage check': 0,  
        'Max extra experiments': 0,
        'Mean steps': 1
    })
    
    # Chinese Remainder method
    method_name = 'Chinese remainder'
    primez = get_chinese_prameters(n_compounds, differentiate)
    methods.append({
        'Pooling strategy': method_name,
        'Max experiments': np.sum(primez),  
        'Max samples per pool': n_compounds/np.min(primez),
        'N pools': np.sum(primez),
        'Percentage check': 0,  
        'Max extra experiments': 0,
        'Mean steps': 1
    })
    
    # Chinese Remainder backtrack method
    method_name = 'Ch. rm. bktrk'
    primez = get_chinese_prameters(n_compounds, differentiate, backtrack=True)
    methods.append({
        'Pooling strategy': method_name,
        'Max experiments': np.sum(primez),  
        'Max samples per pool': n_compounds/np.min(primez),
        'N pools': np.sum(primez),
        'Percentage check': 0,  
        'Max extra experiments': 0,
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
            'Max experiments': t,  
            'Max samples per pool':MS,
            'N pools': t,
            'Percentage check': 0,  
            'Max extra experiments': 0 ,
            'Mean steps':1 
        })
        methods.append({
            'Pooling strategy': method_name,
            'Max experiments': t,  
            'Max samples per pool':MS,
            'N pools': t,
            'Percentage check': 0,  
            'Max extra experiments': 0,
            'Mean steps':1 
        })
        

    method_name = 'Hierarchical'
    metricas=calculate_metrics_hierarchical_fast(n_compounds, differentiate, checks=checks_hierarchical)
    methods.append({
            'Pooling strategy': method_name,
            'Max experiments': metricas[0],  
            'Max samples per pool':metricas[1],
            'N pools': metricas[2],
            'Percentage check': metricas[3],  
            'Max extra experiments': np.round(metricas[4],2),
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


# Utility functions from rand_WA_fast.py

def get_max_C(n_compounds, max_compounds):
    return int(n_compounds/2) if max_compounds==0 else max_compounds

def get_min_C(n_compounds, MC):
    return int(np.sqrt(n_compounds)) if int(np.sqrt(n_compounds))<MC else int(MC/2)

def get_min_W(n_compounds):
    return int(np.log2(n_compounds)*0.5)

def get_max_W(n_compounds):
    return int(2*np.sqrt(n_compounds))

def mean_metrics_fast(well_assigner, differentiate, max_checks=1e0, scaler=1e3, mp=1e-5, **kwargs):
    BT=well_assigner.shape[1]
    n_compounds=well_assigner.shape[0]
    MC=0
    MP=1
    ls_combs=[]
    ls_diffs=[]
    difo=differentiate
    while (MC<max_checks*scaler or MP>mp) and difo>0:
        MCI=math.comb(n_compounds,differentiate)
        MC+=MCI
        mp=MCI/MC
        ls_combs.append(MCI)
        ls_diffs.append(difo)
        difo-=1
        

    counts=[]
    probi=np.array(ls_combs, dtype=float)
    probi/=np.sum(probi)
    differis=np.random.choice(ls_diffs, int(max_checks*scaler), p=probi)
    #differo=differentiate
    #for _ in range(int(max_checks*scaler)):
    MC=np.max([int(max_checks*scaler/10+5), n_compounds+1])
    for differo in differis:
        rnd_pos=np.random.choice(np.arange(n_compounds), differo, replace=False)
        readout=np.any(well_assigner[rnd_pos], axis=0)
        decoded=fast_decode(well_assigner=well_assigner, differentiate=differentiate, 
                            readout=readout, max_checks=MC)
        try:
            decoded_set=set([x for xx in decoded for x in xx])
            NC0=len(decoded_set)
            NC=np.min([NC0,len(decoded)])
        except:
            NC=len(decoded)
        counts.append(NC)
    counts=np.array(counts)
    #ET=np.max(counts-1)#/len(counts)
    ET=np.sum(counts)/len(counts)
    ET =ET if ET>1 else 0#if ET<well_assigner.shape[0] else well_assigner.shape[0]
    ER=np.sum(counts>1)/np.sum(counts>0)
    rounds=ER+1
    p_check=np.round(ER*100)
    return BT+ET, ET,  rounds, p_check
    
 

def decode_precomp(well_assigner:np.array, differentiate:int, 
                   scrambler:dict, readout:np.ndarray, max_differentiate=-1, sweep=False, **kwargs) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    if max_differentiate<1:
        N=well_assigner.shape[0]
        sc_list=np.arange(N).tolist()
        for i in range(differentiate):
            diff=i+1
            if diff ==1:
                full_well_assigner=well_assigner.copy()
            else:
                this_sc=scrambler[diff]
                #print(this_sc)
                #print(well_assigner)
                #print(diff)
                full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
                sc_list.extend(this_sc.tolist())
        #outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=True)
        
        if sweep:
            outcome_dict={}
            outcomes=np.unique(full_well_assigner, axis=0, return_counts=False).astype(int)
            for outcome in outcomes:
                idxs = np.all(outcome == full_well_assigner, axis=1)
                outcome_dict.update({tuple_to_str(tuple(outcome)):list(itertools.compress(sc_list,idxs))})
            return outcome_dict

        else:
            idxs = np.all(readout == full_well_assigner, axis=1)
            return list(itertools.compress(sc_list,idxs))
        
    else:
        full_od={}
        N=well_assigner.shape[0]
        sc_list=np.arange(N).tolist()
        for differentiate in range(max_differentiate):

            diff=differentiate+1
            if diff ==1:
                full_well_assigner=well_assigner.copy()
            else:
                this_sc=scrambler[diff]
                print(this_sc)
                print(well_assigner)
                print(diff)
                full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
                sc_list.extend(this_sc.tolist())
        #outcomes,_=np.unique(full_well_assigner, axis=0, return_counts=True)
        
            if sweep:
                outcome_dict={}
                outcomes=np.unique(full_well_assigner, axis=0, return_counts=False).astype(int)
                for outcome in outcomes:
                    idxs=np.prod(outcome==full_well_assigner, axis=1)
                    outcome_dict.update({tuple_to_str(tuple(outcome)):list(itertools.compress(sc_list,idxs))})
                full_od.update({diff:outcome_dict})

            else:
                idxs=np.prod(readout==full_well_assigner, axis=1)
            full_od.update({diff:list(itertools.compress(sc_list,idxs))})
        return full_od
    
# Random Design Functions from rand_WA_fast.py

def fast_decode(well_assigner:np.array, differentiate:int, readout:np.ndarray, 
                max_checks=1e4, **kwargs):
    WA=well_assigner
    n_pools=WA.shape[1]

    if np.max(readout)>1 or len(readout)!=n_pools:
        readout_bin_ls = [1 if i in readout else 0 for i in range(n_compounds)]
        readout_bl=np.array(readout_bin_ls)
    else:
        readout_bl=readout
    mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
    original_indices = np.where(mask)[0]
    filtered_WA = WA[mask]
    n_compounds=filtered_WA.shape[0]
    if n_compounds<differentiate:
        differentiate=n_compounds
    if n_compounds<2:
        if n_compounds==1:
            decoded=[original_indices[0]] 
        else:
            decoded=[]
    else:
        MC=0
        MP=1
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
        if MC>max_checks:
            decoded = [int(original_indices[idx]) for idx in range(len(original_indices))]
        else:
            scrambler={1:np.arange(n_compounds)}
            for j in range(2,differentiate+1):
                scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
            decoded_pre=decode_precomp(well_assigner=filtered_WA, differentiate= differentiate, scrambler=scrambler, 
                    readout=readout_bl)
            decoders = [combination if isinstance(combination, list) else [combination] for combination in decoded_pre]
            decoded = [[int(original_indices[idx]) for idx in combination] for combination in decoders]
    return decoded

def find_rand_params_precomp(n_compounds:int, n_compounds_per_well=0, n_wells=0, guesses=0, 
                     max_compounds=0, max_redundancy=4, min_redundancy=1, differentiate=1, **kwargs):
    skip_compounds=True
    skip_wells=True
    if n_compounds_per_well==0:
        skip_compounds=False
    if n_wells==0:
        skip_wells=False
    if guesses==0:
        guesses=n_compounds

    MC= get_max_C(n_compounds, max_compounds)
    mc=get_min_C(n_compounds, MC)
    arr_comp=np.arange(int(mc),int(MC+1))
    mw=get_min_W(n_compounds)
    MW=get_max_W(n_compounds)
    while MW-mw<10:
        mw=int(abs(mw-1))
        MW=int(MW+1)

    arr_wells=np.arange(mw,MW+1)
    min_tests=np.inf
    N_tries=0
    for comp in arr_comp:
        if skip_compounds:
            comp=n_compounds_per_well

        for wells in arr_wells:
            if skip_wells:
                if skip_compounds:
                    return 0
                wells=n_wells
                
            if comp*wells>(max_redundancy)*n_compounds*np.log2(n_compounds) or comp*wells<(min_redundancy)*n_compounds*np.log2(n_compounds): continue 
            WA_tmp, mean_exp, p_check=evaluate_rand_design(n_compounds=n_compounds,n_wells=wells,n_compounds_per_well=comp, return_me=True, guesses=guesses, differentiate=differentiate, **kwargs)
            N_tries+=1
            if mean_exp<min_tests:
                Comp=comp
                Wells=wells
                min_tests=mean_exp
                min_wa=WA_tmp
                min_pcheck=p_check
            if skip_wells:
                break
        if skip_compounds:
            break

    print('\n')
    print('----------------------------------------------------------------------------------------------------------')
    print("Evaluated %s different random designs each with %s configurations " % (N_tries, guesses))
    print('\n')

    return Comp, Wells, min_tests, min_wa, min_pcheck

def evaluate_rand_design(n_compounds:int,  differentiate:int, n_compounds_per_well=0, 
                        n_wells=0, guesses=0,  return_me=False, **kwargs):
    min_tests=np.inf
    second_axis=np.tile(np.arange(n_wells),n_compounds_per_well).reshape(n_compounds_per_well,-1)
    for i in range(guesses):
        # Generate idt with no duplicates within each well (column), but duplicates across wells OK
        idt=np.zeros((n_compounds_per_well,n_wells), dtype=int)
        for well_idx in range(n_wells):
            idt[:, well_idx] = np.random.choice(n_compounds, n_compounds_per_well, replace=False)
        well_assigner=np.zeros((n_compounds,n_wells))==1
        well_assigner[idt, second_axis]=True

        mean_exp, _, _, p_check= mean_metrics_fast(well_assigner=well_assigner,
                                                    differentiate=differentiate, **kwargs)
        if p_check<1:
            if return_me:
                return well_assigner,  mean_exp, p_check
            return well_assigner
        elif mean_exp<min_tests: 
            best_wa=well_assigner.copy()
            min_tests=mean_exp
            min_pcheck=p_check

    if return_me:
        return best_wa,  min_tests, min_pcheck
    
    return best_wa

def assign_wells_random_precomp(n_compounds:int,  differentiate:int,n_compounds_per_well=0, 
                        n_wells=0, guesses=0, Evaluate=False, return_me=False, **kwargs)->np.array:

    _,_, min_tests, WA_rand, p_check=find_rand_params_precomp(n_compounds=n_compounds, differentiate=differentiate, 
                                 n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, 
                                 guesses=guesses, **kwargs)
    if return_me:
        return WA_rand,  min_tests, p_check
    
    return WA_rand

def check_Rand_in_WApath(WApath):
    try:
        files = os.listdir(WApath)
        for file in files:
            if file.startswith('WA_Random_N_'):
                return True
        return False
    except Exception as e:
        return str(e)

def append_random_metrics_to_metrics(dpath, N, diff, nw, ms, pc, me, max_dilution):
    """Append random design metrics to the summary CSV file if it exists."""
    summary_file = os.path.join(dpath, f'Metrics_N_{N}_diff_{diff}.csv')
    if not os.path.exists(summary_file):
        return False
    
    try:
        # Read existing CSV
        df = pd.read_csv(summary_file)
        
        if 'Max experiments per sample' not in df.columns:
            df['Max experiments per sample'] = np.nan
        
        # Ensure N and diff columns exist
        if 'N' not in df.columns:
            df['N'] = N
        if 'diff' not in df.columns:
            df['diff'] = diff

        # Check if Random row already exists
        if (df['Pooling strategy'] == 'Random').any():
            # Update existing Random row
            df.loc[df['Pooling strategy'] == 'Random', 'N'] = N
            df.loc[df['Pooling strategy'] == 'Random', 'diff'] = diff
            df.loc[df['Pooling strategy'] == 'Random', 'Mean experiments'] = me
            df.loc[df['Pooling strategy'] == 'Random', 'Max samples per pool'] = ms
            df.loc[df['Pooling strategy'] == 'Random', 'N pools'] = nw
            df.loc[df['Pooling strategy'] == 'Random', 'Percentage check'] = pc
            df.loc[df['Pooling strategy'] == 'Random', 'Max experiments per sample'] = max_dilution
            extra = max(me - nw, 0.0)
            df.loc[df['Pooling strategy'] == 'Random', 'Mean extra experiments'] = extra
            df.loc[df['Pooling strategy'] == 'Random', 'Mean steps'] = 1 + pc/100.0
        else:
            # Append new Random row
            extra = max(me - nw, 0.0)
            random_row = {
                'N': N,
                'diff': diff,
                'Pooling strategy': 'Random',
                'Mean experiments': me,
                'Max samples per pool': ms,
                'N pools': nw,
                'Percentage check': pc,
                'Max experiments per sample': max_dilution,
                'Mean extra experiments': extra,
                'Mean steps': 1 + pc/100.0
            }
            df = pd.concat([df, pd.DataFrame([random_row])], ignore_index=True)
        
        # Round numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        df[numeric_cols] = df[numeric_cols].round(2)
        
        # Reorder columns to have N and diff first
        
        # Write back to CSV
        df.to_csv(summary_file, index=False)
        return True
    except Exception as e:
        print(f"Warning: Could not update summary file {summary_file}: {e}")
        return False

def rand_sweep_diff(n_compounds, max_diff, Npath, max_prev, **kwargs):
    if 'differentiate' in kwargs.keys():
        del kwargs['differentiate']
    
    
    dpath=os.path.join(Npath,'diff_'+str(max_diff))
    WApath=os.path.join(dpath,'WAs')
    if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
        return
    max_diffo= max_diff if max_diff >0 else int(max_prev*n_compounds)+2
    for di in range(max_diffo):
        if di+1>max_prev*n_compounds:
            break
        start_time = time.time()
        diff=di+1
        if diff==1:
            dpath=os.path.join(Npath,'diff_'+str(diff))
            WApath=os.path.join(dpath,'WAs')

            if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
                continue

            WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                            differentiate=diff, return_me=True, **kwargs )
            extra_exp=WA_rand.shape[1]+min_tests

            if kwargs['cleanup']=='one_liner' or kwargs['cleanup']=='full' or kwargs['cleanup']=='True':
                filenames = next(os.walk(WApath), (None, None, []))[2]
                for fname in filenames:
                    if fname.startswith('WA_Random_N_'):
                        os.remove(os.path.join(WApath,fname))

            if kwargs['cleanup']=='WA' or kwargs['cleanup']=='full' or kwargs['cleanup']=='True':
                filenames = next(os.walk(dpath), (None, None, []))[2]
                for fname in filenames:
                    if fname.startswith('Random_diff_'):
                        os.remove(os.path.join(dpath,fname))

            full_file_dir=os.path.join(dpath,'Random_diff_'+str(diff)+'_NS_'+
                                            str(n_compounds)+'_NW_'+str(WA_rand.shape[1])+
                                            '_MS_'+str(np.max(np.sum(WA_rand, axis=0)))+
                                            '_PC_'+ str(int(perc_check)) +'_ME_'+str(np.round(min_tests,2))+".txt")
            if not os.path.exists(dpath):
                os.makedirs(dpath)
            if kwargs['one_liner']:
                open(full_file_dir, 'a').close()

            if not os.path.exists(WApath):
                os.makedirs(WApath)
            thisfile=os.path.join(WApath,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+
                                    '_ME_'+str(np.round(min_tests,2))+'.csv')
            pd.DataFrame(
                WA_rand.astype(int),
                index=[f"Sample_{i}" for i in range(WA_rand.shape[0])],
                columns=[f"Pool_{j}" for j in range(WA_rand.shape[1])],
            ).to_csv(thisfile)

            # Append to summary CSV if it exists
            append_random_metrics_to_metrics(
                dpath=dpath,
                N=n_compounds,
                diff=diff,
                nw=WA_rand.shape[1],
                ms=np.max(np.sum(WA_rand, axis=0)),
                pc=int(perc_check),
                me=np.round(min_tests, 2),
                max_dilution=int(np.max(np.sum(WA_rand, axis=1)))
            )

        else:
            dpath=os.path.join(Npath,'diff_'+str(diff))
            WApath=os.path.join(dpath,'WAs')
            
            if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
                continue

            WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                            differentiate=diff, return_me=True, **kwargs  )
            extra_exp=WA_rand.shape[1]+min_tests

            if kwargs['cleanup']=='one_liner' or kwargs['cleanup']=='full' or kwargs['cleanup']=='True':
                filenames = next(os.walk(WApath), (None, None, []))[2]
                for fname in filenames:
                    if fname.startswith('WA_Random_N_'):
                        os.remove(os.path.join(WApath,fname))

            if kwargs['cleanup']=='WA' or kwargs['cleanup']=='full' or kwargs['cleanup']=='True':
                filenames = next(os.walk(dpath), (None, None, []))[2]
                for fname in filenames:
                    if fname.startswith('Random_diff_'):
                        os.remove(os.path.join(dpath,fname))

            full_file_dir=os.path.join(dpath,'Random_diff_'+str(diff)+'_NS_'+
                                            str(n_compounds)+'_NW_'+str(WA_rand.shape[1])+
                                            '_MS_'+str(np.max(np.sum(WA_rand, axis=0)))+
                                            '_PC_'+ str(int(perc_check)) +'_ME_'+str(np.round(min_tests,2))+".txt")

            if not os.path.exists(dpath):
                os.makedirs(dpath)

            if kwargs['one_liner']:
                open(full_file_dir, 'a').close()
            if not os.path.exists(WApath):
                os.makedirs(WApath)
            thisfile=os.path.join(WApath,'WA_Random_N_'+str(n_compounds)+'_diff_'+str(diff)+
                                    '_ME_'+str(np.round(min_tests,2))+'.csv')
            pd.DataFrame(
                WA_rand.astype(int),
                index=[f"Sample_{i}" for i in range(WA_rand.shape[0])],
                columns=[f"Pool_{j}" for j in range(WA_rand.shape[1])],
            ).to_csv(thisfile)

            # Append to summary CSV if it exists
            append_random_metrics_to_metrics(
                dpath=dpath,
                N=n_compounds,
                diff=diff,
                nw=WA_rand.shape[1],
                ms=np.max(np.sum(WA_rand, axis=0)),
                pc=int(perc_check),
                me=np.round(min_tests, 2),
                max_dilution=int(np.max(np.sum(WA_rand, axis=1)))
            )

        DTS=np.round((time.time() - start_time),2)
        DTD=DTS//86400
        DTH=DTS//3600-DTD*24
        DTM=DTS//60-DTH*60-DTD*24*60
        DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
        print("%s days %s hours %s minutes and %s seconds required for N= %s and differentiate %s" % 
                (DTD, DTH, DTM, DTS, n_compounds, diff))
        print('----------------------------------------------------------------------------------------------------------') 


def rand_N_sweep(start, stop, step,dir_WAs, **kwargs):
    n_compounds=start
    while n_compounds<=stop:
        start_time = time.time()
        Npath=os.path.join(dir_WAs,'N_'+str(n_compounds))
        rand_sweep_diff(n_compounds=n_compounds, Npath=Npath, **kwargs)
        DTS=np.round((time.time() - start_time),2)
        DTD=DTS//86400
        DTH=DTS//3600-DTD*24
        DTM=DTS//60-DTH*60-DTD*24*60
        DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
        print('\n')
        print('\n')
        print('##########################################################################################################')
        print('##########################################################################################################')
        print('\n')
        print("%s days %s hours %s minutes and %s seconds overall required for N= %s" % (DTD, DTH, DTM, DTS, n_compounds))
        print('\n')
        print('##########################################################################################################')
        print('##########################################################################################################')

        n_compounds+=step


def get_wa_filename(save_dir, n_compounds, diff, method):
    """Generate consistent filename pattern"""
    return os.path.join(
        save_dir,
        f'N_{n_compounds}',
        f'diff_{diff}',
        'WAs',
        f'WA_{method}_N_{n_compounds}_diff_{diff}.csv'
    )

def process_n_compounds(**kwargs):
    """
    Processes WA computations for a specific n_compounds value with diff handling
    and file existence checks
    """

    n_compounds = kwargs['n_compounds']
    max_diff = kwargs['max_diff']
    save_dir = kwargs['save_dir']
    max_dims = kwargs['max_dims']
    max_prev = kwargs['max_prev']
    timeit = kwargs.get('timeit', False)
    
    # Cleanup non-random WA files if requested - do this at the very beginning
    cleanup = kwargs.get('cleanup', False)
    if cleanup in ['True', 'true', True, 'WA', 'full', 'one_liner']:
        n_dir = os.path.join(save_dir, f'N_{n_compounds}')
        if os.path.exists(n_dir):
            for diff_entry in os.listdir(n_dir):
                if not diff_entry.startswith('diff_'):
                    continue
                diff_dir = os.path.join(n_dir, diff_entry, 'WAs')
                if os.path.exists(diff_dir):
                    for fname in os.listdir(diff_dir):
                        if fname.endswith('.csv') and not fname.startswith('WA_Random_N_'):
                            os.remove(os.path.join(diff_dir, fname))
                            if timeit:
                                print(f"Removed {fname} from {diff_dir}")
    
    # Generate multidim methods dynamically
    multidim_methods = []
    for i in np.arange(3, int(np.ceil(np.log(n_compounds)))+1):
        if i > max_dims:
            continue
        multidim_methods.append(f'multidim-{i}')
    
    # Compute diff-independent methods only if needed
    WA_list = []
    methods = []
    computed_diff_independent = False
    
    # Check if all diff-independent files exist for diff=1
    all_exist = True
    for method in multidim_methods + ['Matrix', 'Binary']:
        if not os.path.exists(get_wa_filename(save_dir, n_compounds, 1, method)):
            all_exist = False
            break
    
    # Compute diff-independent methods if any are missing
    if not all_exist:
        computed_diff_independent = True
        # 1. Compute multidim methods
        for method in multidim_methods:
            dim = int(method.split('-')[1])
            WA = assign_wells_multidim(n_dims=dim, **kwargs)
            WA_list.append(WA)
            methods.append(method)
        
        # 2. Compute Matrix and Binary
        WA_mat = assign_wells_mat(**kwargs)
        WA_list.append(WA_mat)
        methods.append('Matrix')
        
        WA_bin = assign_wells_bin(**kwargs)
        WA_list.append(WA_bin)
        methods.append('Binary')
    elif timeit:
        print(f"All diff-independent files exist for n={n_compounds}, skipping computation")
    
    # Process each differentiation level
    max_diffo= max_diff if max_diff >0 else int(max_prev*n_compounds)+1
    for diff in range(1, max_diffo + 1):
        if diff>max_prev*n_compounds:
            break
        current_kwargs = kwargs.copy()
        current_kwargs['differentiate'] = diff
        
        # Create directory structure
        diff_dir = os.path.join(save_dir, f'N_{n_compounds}', f'diff_{diff}', 'WAs')
        os.makedirs(diff_dir, exist_ok=True)
        
        # Handle diff-independent methods
        for method in multidim_methods + ['Matrix', 'Binary']:
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                if diff == 1 and computed_diff_independent:
                    # Find method index and save
                    idx = methods.index(method)
                    pd.DataFrame(
                        WA_list[idx].astype(int),
                        index=[f"Sample_{i}" for i in range(WA_list[idx].shape[0])],
                        columns=[f"Pool_{j}" for j in range(WA_list[idx].shape[1])],
                    ).to_csv(dst_file)
                    if timeit:
                        print(f"Saved {method} for diff={diff}")
                else:
                    # Try to copy from diff=1
                    src_file = get_wa_filename(save_dir, n_compounds, 1, method)
                    if os.path.exists(src_file):
                        shutil.copy(src_file, dst_file)
                        if timeit:
                            print(f"Copied {method} to diff={diff}")
                    elif timeit:
                        print(f"Warning: Base file missing for {method} at diff={diff}")
        
        # Handle diff-dependent methods
        for method in ['STD', 'Chinese remainder', 'Ch. rm. bktrk']:
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                # Compute only if file doesn't exist
                if method == 'STD':
                    WA = assign_wells_STD(**current_kwargs)
                elif method == 'Chinese remainder':
                    WA = assign_wells_chinese(**current_kwargs)
                else:  # Ch. rm. bktrk
                    WA = assign_wells_chinese(backtrack=True, **current_kwargs)
                
                    pd.DataFrame(
                        WA.astype(int),
                        index=[f"Sample_{i}" for i in range(WA.shape[0])],
                        columns=[f"Pool_{j}" for j in range(WA.shape[1])],
                    ).to_csv(dst_file)
                if timeit:
                    print(f"Computed {method} for diff={diff}")
            elif timeit:
                print(f"Skipping {method} for diff={diff} (already exists)")
        
        # Handle Chinese special (only for diff 2 or 3)
        if diff in [2, 3]:
            method = 'Chinese special'
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                WA = assign_wells_chinese(special_diff=True, **current_kwargs)
                pd.DataFrame(
                    WA.astype(int),
                    index=[f"Sample_{i}" for i in range(WA.shape[0])],
                    columns=[f"Pool_{j}" for j in range(WA.shape[1])],
                ).to_csv(dst_file)
                if timeit:
                    print(f"Computed {method} for diff={diff}")
            elif timeit:
                print(f"Skipping {method} for diff={diff} (already exists)")

def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    """
    Main loop to process all n_compounds values
    """
    current = start
    while current <= stop:
        if kwargs.get('timeit'):
            print(f"Processing {current} compounds")
            time0 = time.time()
        
        # Process this n_compounds value
        kwargs['n_compounds'] = current
        process_n_compounds(**kwargs)
        
        if kwargs.get('timeit'):
            DTS=np.round((time.time() - time0),2)
            DTD=DTS//86400
            DTH=DTS//3600-DTD*24
            DTM=DTS//60-DTH*60-DTD*24*60
            DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
            print('\n')
            print("%s days %s hours %s minutes and %s seconds required for N= %s compounds" % 
                  (DTD, DTH, DTM, DTS, current))
            print('\n')

        current += step


# Function to generate expected error table for prevalence tab
def expected_error_table(N, prevalence, max_diff=4, min_N=5, correct=False, max_error=0.05, extra_steps=2):
 # or set to N or another reasonable value
    diff_values = np.arange(1, int(max_diff)+extra_steps+1)
    extreme_diff= int(np.ceil(scipy.stats.binom.isf(max_error, N, prevalence)))

    if extreme_diff>max_diff+extra_steps:
        rooty=np.power(extreme_diff/max_diff,1/(extra_steps))
        for i in range(extra_steps):
            new_diff=int(np.ceil(max_diff*np.power(rooty, i+1)))
            diff_values[max_diff+i]=new_diff
    #elif max_diff<extreme_diff<max_diff+extra_steps:
    #    diff_values = np.arange(1, int(extreme_diff)+1)
    else:
        diff_values = np.arange(1, int(np.max([max_diff, extreme_diff]))+1)

    pop_sizes = []
    N_current = int(np.ceil(N))
    factor=[]
    N_factor=[]
    NF=1
    # Build the population sizes as described, always take the ceiling to be on the safe side
    while N_current >= min_N:
        pop_sizes.append(N_current)
        factor.append(NF)
        N_factor.append(str(NF)+' X '+str(N_current))
        N_current = int(np.ceil(N_current / 2))
        NF=int(NF*2)
        if N_current < min_N:
            break
        N_factor.append(str(NF)+' X '+str(N_current))
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2.5)
        N_current = int(np.ceil(N_current / 2.5))
        if N_current < min_N:
            break
        N_factor.append(str(NF)+' X '+str(N_current))
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2)
        N_current = int(np.ceil(N_current / 2))
    # Remove duplicates and sort descending
    #pop_sizes = sorted(set([x for x in pop_sizes if x >= min_N]), reverse=True)
    
    # Build the table
    
    data = []
    CF=pop_sizes
    for N_val,NF in zip(pop_sizes,factor):
        row = []
        for diff in diff_values:
            p_error = scipy.stats.binom.sf(diff, N_val, prevalence)
            if correct:
                p_error=1-(1-p_error)**int(NF)
                CF=N_factor
            if p_error<1e-15:
                p_error=0
            row.append(p_error)
        data.append(row)
    df = pd.DataFrame(data, index=CF, columns=diff_values)
    # No minimum error logic needed
    df.index.name = "N"
    df.columns.name = "Differentiate"
    return df





def parse_method(filename: str) -> str:
	match = re.match(r"^WA_(?P<method>.+)_N_\d+_diff_\d+.*\.csv$", filename)
	return match.group('method') if match else 'unknown'


def iter_wa_files(base_dir: str):
	for n_entry in sorted(os.listdir(base_dir)):
		if not n_entry.startswith('N_'):
			continue
		n_dir = os.path.join(base_dir, n_entry)
		if not os.path.isdir(n_dir):
			continue
		try:
			n_val = int(n_entry.split('_', 1)[1])
		except ValueError:
			continue

		for diff_entry in sorted(os.listdir(n_dir)):
			if not diff_entry.startswith('diff_'):
				continue
			diff_dir = os.path.join(n_dir, diff_entry)
			if not os.path.isdir(diff_dir):
				continue
			try:
				diff_val = int(diff_entry.split('_', 1)[1])
			except ValueError:
				continue

			wa_dir = os.path.join(diff_dir, 'WAs')
			if not os.path.isdir(wa_dir):
				continue

			for fname in sorted(os.listdir(wa_dir)):
				if not fname.endswith('.csv'):
					continue
				yield n_val, diff_val, parse_method(fname), os.path.join(wa_dir, fname)


def compute_metrics(n_compounds: int, diff: int, wa: np.ndarray):
	mean_exp, extra, rounds, p_check = mean_metrics_fast(well_assigner=wa, differentiate=diff)
	nw = wa.shape[1]
	ms = int(np.max(np.sum(wa, axis=0)))
	max_dilution = int(np.max(np.sum(wa, axis=1)))+1
	me = float(np.round(mean_exp, 2))
	pc = int(p_check)
	return {
		'N': n_compounds,
		'diff': diff,
		'Mean experiments': me,
		'Max samples per pool': ms,
		'N pools': nw,
		'Percentage check': pc,
		'Max experiments per sample': max_dilution,
		'Mean extra experiments': float(np.round(extra, 2)),
		'Mean steps': rounds
	}


def series_to_readout_list(ser: pd.Series, n_pools_local: int):
    vals = [str(v).strip() for v in ser.tolist() if pd.notna(v) and str(v).strip() != ""]
    if len(vals) == 0:
        return None, "Empty row"
    # Case: single cell input
    if len(vals) == 1:
        cell = vals[0]
        # comma-separated
        if "," in cell:
            parts = [p.strip() for p in cell.split(",") if p.strip() != ""]
            if any(not p.isdigit() for p in parts):
                return None, f"Invalid token(s) in list: {parts}"
            return [int(p) for p in parts], None
        # plain 0/1 string
        if set(cell).issubset({"0", "1"}) and len(cell) == n_pools_local:
            return [int(ch) for ch in cell], None
        # space-separated integers
        parts = [p for p in re.split(r"[\s;]+", cell) if p]
        if parts and all(p.isdigit() for p in parts):
            return [int(p) for p in parts], None
        return None, f"Unrecognized single-cell format: '{cell}'"
    # Case: multi-column row
    try:
        ints = [int(float(v)) for v in vals]
    except Exception:
        return None, f"Non-integer values across columns: {vals}"
    uniq = set(ints)
    if uniq.issubset({0, 1}) and len(ints) == n_pools_local:
        return ints, None
    return ints, None


def normalize_readout_df(readout_df: pd.DataFrame) -> pd.DataFrame:
    if readout_df is None:
        return pd.DataFrame()
    out = readout_df.copy()
    out = out.map(lambda x: x.strip() if isinstance(x, str) else x)
    out = out.replace(r"^\s*$", np.nan, regex=True)
    out = out.dropna(axis=0, how='all').dropna(axis=1, how='all')
    return out


def is_int_like_value(val) -> bool:
    try:
        if pd.isna(val):
            return False
        fv = float(str(val).strip())
        if not np.isfinite(fv):
            return False
        return abs(fv - round(fv)) < 1e-9
    except Exception:
        return False


def _tokens_from_payload(payload):
    if isinstance(payload, pd.Series):
        raw_vals = payload.tolist()
    elif isinstance(payload, (list, tuple, np.ndarray)):
        raw_vals = list(payload)
    else:
        raw_vals = [payload]

    vals = []
    for v in raw_vals:
        if pd.isna(v):
            continue
        sv = str(v).strip()
        if sv == "":
            continue
        vals.append(sv)

    if len(vals) == 1:
        cell = vals[0]
        if re.search(r"[\s,;]", cell):
            parts = [p.strip() for p in re.split(r"[\s,;]+", cell) if p.strip() != ""]
            return parts
    return vals


def classify_single_readout_payload(payload, n_pools_local: int):
    tokens = _tokens_from_payload(payload)
    if len(tokens) == 0:
        return None, None, "Error: empty readout."

    try:
        nums = [float(t) for t in tokens]
    except Exception:
        # Try to parse as a single comma-separated string of numerics
        if len(tokens) == 1:
            csv_str = str(tokens[0]).strip()
            try:
                nums = [float(x.strip()) for x in csv_str.split(',') if x.strip()]
                tokens = [str(x) for x in nums]
                if len(nums) == 0:
                    return None, None, "Error: empty readout."
            except Exception:
                return None, None, "Error: readout must be numeric for continuous mode, or integers for binary mode."
        else:
            return None, None, "Error: readout must be numeric for continuous mode, or integers for binary mode."

    arr_float = np.array(nums, dtype=float)
    int_like = [is_int_like_value(v) for v in nums]

    # Condition 1: length == n_pools -> continuous decode (unless it's a strict binary vector)
    if len(nums) == n_pools_local:
        if all(int_like):
            ints = [int(round(v)) for v in nums]
            # Only treat as binary if it's a 0/1 binary vector
            if set(ints).issubset({0, 1}):
                return "binary", np.array(ints, dtype=int), None
        # Otherwise it's continuous (even if values are integers)
        return "continuous", arr_float, None

    # Condition 2 (alternative): integer list with fewer values and max <= n_pools -> binary decode.
    # This catches pool indices like [0, 2, 4] or [1, 3]
    if all(int_like) and len(nums) < n_pools_local:
        ints = [int(round(v)) for v in nums]
        if min(ints) >= 0 and max(ints) <= n_pools_local:
            return "binary", np.array(sorted(ints), dtype=int), None

    return None, None, (
        f"Error: input does not satisfy either condition. "
        f"Expected either {n_pools_local} numeric values (continuous/non-binary or binary), "
        "or an integer list with max <= number of pools."
    )


def decode_single_readout_payload(payload, n_pools_local: int, WA: np.ndarray, diff: int, readout_id: str = None, min_signal=None, diluting=True, alpha=1.0, l1_ratio=0.5):
    mode, arr, err = classify_single_readout_payload(payload, n_pools_local)
    if err:
        out = {
            "decoded_type": "error",
            "decoder_output": err,
        }
        if readout_id is not None:
            out["readout_id"] = readout_id
        return out

    if mode == "continuous":
        decoded_type, decoder_output= decode_continuous_lasso(
            arr,
            WA,
            diff=diff,
            min_signal=min_signal,
            diluting=diluting,
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        out = {
            "decoded_type": decoded_type,
            "decoder_output": decoder_output,
        }
        if readout_id is not None:
            out["readout_id"] = readout_id
        return out

    decoded_type, decoder_output, _, _, _, _ = decode_with_filter(arr, WA, diff)
    out = {
        "decoded_type": decoded_type,
        "decoder_output": decoder_output,
    }
    if readout_id is not None:
        out["readout_id"] = readout_id
    return out


def _row_decode_binary_like(values, n_pools_local: int, WA: np.ndarray, diff: int):
    ints = []
    for v in values:
        if pd.isna(v):
            continue
        sv = str(v).strip()
        if sv == "":
            continue
        if not is_int_like_value(sv):
            return {
                "decoded_type": "error",
                "decoder_output": f"Non-integer value found in binary/integer readout: '{sv}'",
            }
        iv = int(round(float(sv)))
        if iv < 0 or iv > n_pools_local:
            return {
                "decoded_type": "error",
                "decoder_output": f"Entries exceed number of pools ({n_pools_local}): {[iv]}",
            }
        ints.append(iv)

    if len(ints) == 0:
        return {
            "decoded_type": "error",
            "decoder_output": "Error: empty readout row.",
        }

    # If row is exactly a binary vector, keep as vector; otherwise treat as pool-index list.
    if len(ints) == n_pools_local and set(ints).issubset({0, 1}):
        readout_arr = np.array(ints, dtype=int)
    else:
        readout_arr = np.array(sorted(ints), dtype=int)

    decoded_type, decoder_output, _, _, _, _ = decode_with_filter(readout_arr, WA, diff)
    return {
        "decoded_type": decoded_type,
        "decoder_output": decoder_output,
    }


def _is_condition1_matrix(readout_df: pd.DataFrame, n_pools_local: int) -> bool:
    if readout_df is None or readout_df.empty:
        return False
    for _, row in readout_df.iterrows():
        tokens = _tokens_from_payload(row)
        if len(tokens) != n_pools_local:
            return False
        try:
            nums = [float(t) for t in tokens]
        except Exception:
            return False
        if all(is_int_like_value(v) for v in nums):
            ints = [int(round(v)) for v in nums]
            if set(ints).issubset({0, 1}):
                return False
    return True


def _is_condition2_binary_vector_matrix(readout_df: pd.DataFrame, n_pools_local: int) -> bool:
    if readout_df is None or readout_df.empty:
        return False
    for _, row in readout_df.iterrows():
        tokens = _tokens_from_payload(row)
        if len(tokens) != n_pools_local:
            return False
        try:
            nums = [float(t) for t in tokens]
        except Exception:
            return False
        if not all(is_int_like_value(v) for v in nums):
            return False
        ints = [int(round(v)) for v in nums]
        if not set(ints).issubset({0, 1}):
            return False
    return True


def is_binary_or_integer_matrix(readout_df: pd.DataFrame, n_pools_local: int) -> bool:
    if readout_df is None or readout_df.empty:
        return False
    found_any = False
    for v in readout_df.to_numpy().flatten():
        if pd.isna(v):
            continue
        sv = str(v).strip()
        if sv == "":
            continue
        if not is_int_like_value(sv):
            return False
        iv = int(round(float(sv)))
        if iv < 0 or iv > n_pools_local:
            return False
        found_any = True
    return found_any


def is_continuous_with_id_matrix(readout_df: pd.DataFrame, n_pools_local: int) -> bool:
    if readout_df is None or readout_df.empty:
        return False
    if readout_df.shape[1] != n_pools_local + 1:
        return False

    for _, row in readout_df.iterrows():
        row_id = row.iloc[0]
        if pd.isna(row_id) or str(row_id).strip() == "":
            return False
        vals = pd.to_numeric(row.iloc[1:], errors='coerce')
        if vals.isna().any():
            return False
    return True


def _decode_fallback_row(row: pd.Series, row_idx: int, n_pools_local: int, WA: np.ndarray, diff: int, min_signal=None, diluting=True, alpha=1.0, l1_ratio=0.5):
    tokens = [v for v in row.tolist() if not pd.isna(v) and str(v).strip() != ""]
    if len(tokens) == 0:
        return {
            "readout_id": f"Readout {row_idx + 1}",
            "decoded_type": "error",
            "decoder_output": "Error: empty row.",
        }

    # If first token is non-numeric and there is payload after it, treat as row id.
    first = str(tokens[0]).strip()
    has_id = False
    try:
        float(first)
    except Exception:
        has_id = len(tokens) > 1

    if has_id:
        row_id = first
        payload = tokens[1:]
    else:
        row_id = f"Readout {row_idx + 1}"
        payload = tokens

    return decode_single_readout_payload(
        payload,
        n_pools_local,
        WA,
        diff,
        readout_id=row_id,
        min_signal=min_signal,
        diluting=diluting,
        alpha=alpha,
        l1_ratio=l1_ratio,
    )


def decode_multi_readout_df(readout_df: pd.DataFrame, WA: np.ndarray, diff: int, min_signal=None, diluting=True, readout_ids=None, alpha=1.0, l1_ratio=0.5) -> pd.DataFrame:
    n_pools_local = WA.shape[1]
    base_df = normalize_readout_df(readout_df)
    if base_df is None or base_df.empty:
        return pd.DataFrame([{
            "readout_id": "Readout 1",
            "decoded_type": "error",
            "decoder_output": "Error: Readout CSV appears to be empty.",
        }])

    # Try full-matrix checks on a small set of candidates to handle optional header/index rows.
    candidates = [("raw", base_df)]
    if base_df.shape[0] > 1:
        candidates.append(("drop_first_row", base_df.iloc[1:, :].reset_index(drop=True)))
    if base_df.shape[1] > 1:
        candidates.append(("drop_first_col", base_df.iloc[:, 1:].reset_index(drop=True)))
    if base_df.shape[0] > 1 and base_df.shape[1] > 1:
        candidates.append(("drop_first_row_col", base_df.iloc[1:, 1:].reset_index(drop=True)))

    mode = None
    chosen_df = None

    for _, cand in candidates:
        cand_norm = normalize_readout_df(cand)
        if cand_norm.empty:
            continue
        # Requested precedence:
        # 1) all rows length == n_pools and not binary -> continuous
        if _is_condition1_matrix(cand_norm, n_pools_local):
            mode = "continuous"
            chosen_df = cand_norm
            break
        # 2) all rows length == n_pools and binary OR integer matrix with max <= n_pools -> binary
        if _is_condition2_binary_vector_matrix(cand_norm, n_pools_local) or is_binary_or_integer_matrix(cand_norm, n_pools_local):
            mode = "binary"
            chosen_df = cand_norm
            break

    rows = []
    if mode == "continuous" and chosen_df is not None:
        for i, (_, row) in enumerate(chosen_df.iterrows()):
            out = decode_single_readout_payload(
                row.tolist(),
                n_pools_local,
                WA,
                diff,
                min_signal=min_signal,
                diluting=diluting,
                alpha=alpha,
                l1_ratio=l1_ratio,
            )
            # Use provided readout_ids if available, otherwise use default naming
            if readout_ids and i < len(readout_ids):
                out["readout_id"] = readout_ids[i]
            else:
                out["readout_id"] = f"Readout {i + 1}"
            rows.append(out)
        return pd.DataFrame(rows)

    if mode == "binary" and chosen_df is not None:
        for i, (_, row) in enumerate(chosen_df.iterrows()):
            out = decode_single_readout_payload(
                row.tolist(),
                n_pools_local,
                WA,
                diff,
                min_signal=min_signal,
                diluting=diluting,
                alpha=alpha,
                l1_ratio=l1_ratio,
            )
            # Use provided readout_ids if available, otherwise use default naming
            if readout_ids and i < len(readout_ids):
                out["readout_id"] = readout_ids[i]
            else:
                out["readout_id"] = f"Readout {i + 1}"
            rows.append(out)
        return pd.DataFrame(rows)

    # Fallback: parse row-by-row as user-provided string with id.
    fallback_df = base_df
    for i, (_, row) in enumerate(fallback_df.iterrows()):
        result = _decode_fallback_row(
            row,
            i,
            n_pools_local,
            WA,
            diff,
            min_signal=min_signal,
            diluting=diluting,
            alpha=alpha,
            l1_ratio=l1_ratio,
        )
        # Use provided readout_ids if available
        if readout_ids and i < len(readout_ids):
            result["readout_id"] = readout_ids[i]
        rows.append(result)

    return pd.DataFrame(rows)


def decode_continuous(readout_arr, WA, diff=None, min_signal=None, diluting=True, alpha=0.05, l1_ratio=1):
    try:
        n_compounds, n_pools = WA.shape
        readout_arr = np.array(readout_arr, dtype=float).flatten()
        
        if len(readout_arr) != n_pools:
            return "error", f"Readout length {len(readout_arr)} doesn't match pool count {n_pools}", None, None, n_compounds, None
        
        # Feature matrix: X = WA.T (n_pools, n_compounds)
        X = WA.T.astype(float)
        y = readout_arr
        
        # Objective function: elastic net loss with MSE
        def objective(coef):
            # Predictions
            y_pred = X @ coef
            # MSE loss
            mse_loss = np.sum((y - y_pred) ** 2) / len(y)
            # L1 penalty (lasso part)
            l1_penalty = np.sum(np.abs(coef))
            # L2 penalty (ridge part)
            l2_penalty = np.sum(coef ** 2)
            # Elastic net: alpha * (l1_ratio * L1 + (1 - l1_ratio) * L2)
            elastic_penalty = alpha * (l1_ratio * l1_penalty + (1 - l1_ratio) * l2_penalty)
            
            return mse_loss + elastic_penalty
        
        # Bounds: all coefficients >= 0
        bounds = [(0, None) for _ in range(n_compounds)]
        
        # Initial guess using NNLS (non-negative least squares)
        coefficients_init, _ = nnls(X, y)
        
        # Optimize with non-negativity constraints
        result = minimize(
            objective,
            coefficients_init,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 10000, 'ftol': 1e-6}
        )
        
        # Get the optimized coefficients
        if result.success:
            coefficients = result.x
        else:
            # Fallback to NNLS if optimization fails
            coefficients, _ = nnls(X, y)
        
        # Final enforcement: ensure all are non-negative (handle numerical errors)
        coefficients = np.maximum(coefficients, 0.0)
        
        # Double-check: ensure all coefficients are truly non-negative
        assert np.all(coefficients >= -1e-12), f"Non-negativity constraint violated! Min value: {np.min(coefficients)}"
        
        # Calculate fit quality metrics
        y_pred = X @ coefficients
        residuals = y - y_pred
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2_score = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        # Return the compound concentration estimates
        decoder_output = coefficients.tolist()
        
        fit_info = {
            'r2_score': float(r2_score),
            'coefficients': decoder_output,
            'predictions': y_pred.tolist(),
            'residuals': residuals.tolist(),
            'alpha': alpha,
            'l1_ratio': l1_ratio,
            'n_nonzero': int(np.sum(coefficients > 1e-6)),
            'all_nonnegative': bool(np.all(coefficients >= -1e-10))
        }
        
        return "continuous", decoder_output, coefficients, r2_score, n_compounds, fit_info
    
    except Exception as e:
        return "error", f"Error in decode_continuous: {str(e)}", None, None, None, None





def process_readout_row(row, n_pools, WA, diff, min_signal=None, diluting=True):
    # Check if row is continuous (array of numeric values, not binary/discrete)
    if len(row) == n_pools and set(row) != {0, 1} and set(row) != {True, False}:
        try:
            readout_arr = np.array(row, dtype=float)
            decoded_type, decoder_output = decode_continuous_lasso(
                readout_arr,
                WA,
                diff=diff,
                min_signal=min_signal,
                diluting=diluting,
            )
            return {
                "decoded_type": decoded_type,
                "decoder_output": decoder_output,
            }
        except Exception as e:
            return {
                "decoded_type": "error",
                "decoder_output": f"Error processing continuous readout: {str(e)}",
            }
    
    # Otherwise, process as discrete readout
    rl, err = series_to_readout_list(row, n_pools)
    if rl is None:
        return {
            "decoded_type": "error",
            "decoder_output": f"Error: {err}",
        }
    
    try:
        rl_sorted = sorted([int(x) for x in rl])
    except Exception:
        return {
            "decoded_type": "error",
            "decoder_output": "Error: could not coerce values to integers",
        }
    
    readout_arr = np.array(rl_sorted, dtype=int)
    if len(readout_arr) == 0:
        return {
            "decoded_type": "error",
            "decoder_output": "Error: empty readout after parsing",
        }
    
    if readout_arr.max() > n_pools:
        invalid_entries = [int(val) for val in rl_sorted if int(val) > n_pools]
        msg = (
            f"Entries exceed number of pools ({n_pools}): {invalid_entries}. "
            "Please correct your readout."
        )
        return {
            "decoded_type": "error",
            "decoder_output": msg,
        }
    
    decoded_type, decoder_output, _, _, _, _ = decode_with_filter(readout_arr, WA, diff)
    return {
        "decoded_type": decoded_type,
        "decoder_output": decoder_output,
    }


############## Decoder specific function


def html_to_text(html_str: str) -> str:
    """Convert HTML to plain text by removing HTML tags."""
    text = re.sub(r'<[^>]+>', '', html_str)
    text = text.replace('<br>', '\n').replace('<br/>', '\n').replace('<br />', '\n')
    return text


def cell_is_valid_token(val, n_pools_local: int):
    try:
        # Skip NaN/None/empty
        if pd.isna(val):
            return None
        s = str(val).strip()
        if s == "" or s.lower() == "nan":
            return None
        # integer scalar
        try:
            iv = int(float(s))
            return 0 <= iv <= n_pools_local
        except Exception:
            pass
        # binary string of length n_pools
        if set(s).issubset({"0", "1"}) and len(s) == n_pools_local:
            return True
        # delimited list of ints
        parts = [p for p in re.split(r"[\s,;]+", s) if p]
        if parts and all(p.isdigit() and 0 <= int(p) <= n_pools_local for p in parts):
            return True
        return False
    except Exception:
        return False


def parse_text_readout(readout_str: str):
    if not readout_str or not readout_str.strip():
        return None, "Error: Please enter a readout as comma-separated text."
    readout_list = []
    try:
        for x in readout_str.split(','):
            x_clean = x.strip()
            if x_clean == '':
                continue
            if not x_clean.isdigit():
                return None, f"Error: Invalid entry '{x_clean}' in readout. Please enter only integers separated by commas."
            readout_list.append(int(x_clean))
    except Exception:
        return None, "Error: Could not parse readout. Please enter only integers separated by commas."
    if len(readout_list) == 0:
        return None, "Error: Readout appears to be empty."
    return readout_list, None


def decode_with_filter(readout_arr: np.ndarray, WA: np.ndarray, diff_deco: int):
    n_pools_local = WA.shape[1]
    n_compounds_local = WA.shape[0]

    if np.max(readout_arr) > 1 or len(readout_arr) != n_pools_local:
        readout_bin_ls = [1 if i in readout_arr else 0 for i in range(n_pools_local)]
        readout_use = np.array(readout_bin_ls)
    else:
        readout_use = readout_arr

    readout_bl = np.array(readout_use.astype(bool).astype(int))
    mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
    original_indices = np.where(mask)[0]
    filtered_WA = WA[mask]
    n_compounds_f = filtered_WA.shape[0]
    dval = min(diff_deco, n_compounds_f)

    if n_compounds_f < 2:
        if n_compounds_f == 1:
            decoded = [int(original_indices[0])]
            return "unique", decoded, decoded, None, n_compounds_f, None
        return "error", "No matches found. Check input or increase the differentiate value.", [], None, n_compounds_f, None

    ls_combs = [math.comb(n_compounds_f, i) for i in range(dval)]
    max_combs = np.sum(ls_combs)
    if max_combs > 1e4:
        decoded = [int(original_indices[j]) for j in range(len(original_indices))]
        decoded = sorted(decoded)
        return "putative", decoded, decoded, decoded, n_compounds_f, "max_combs"

    scrambler = {1: np.arange(n_compounds_f)}
    for j in range(2, dval + 1):
        scrambler[j] = np.array(list(itertools.combinations(np.arange(n_compounds_f), j)))

    decoded_pre = decode_precomp(
        well_assigner=filtered_WA,
        differentiate=dval,
        scrambler=scrambler,
        readout=readout_bl,
    )
    decoders = [list(c) if isinstance(c, (list, tuple, np.ndarray)) else [c] for c in decoded_pre]
    decoded = [[int(original_indices[k]) for k in comb] for comb in decoders]

    if len(decoded) == 0:
        return "error", "No matches found. Check input or increase the differentiate value.", decoded, None, n_compounds_f, None
    if len(decoded) == 1:
        return "unique", decoded[0], decoded, None, n_compounds_f, None
    if len(decoded) > n_compounds_f:
        decoded_set = sorted(list(set([x for combo in decoded for x in combo])))
        return "putative", decoded_set, decoded, decoded_set, n_compounds_f, "too_many_combos"
    return "multiple", decoded, decoded, None, n_compounds_f, None


def build_text_message(file_name: str, readout_list: list, diff_deco: int, n_compounds: int, n_pools: int,
                        decoded_type: str, decoded: list, decoded_set: list, n_compounds_f: int, putative_reason: str):
    lines = []
    lines.append(f"Processing file {file_name} with max. {diff_deco} positive samples and readout:")
    lines.append(f"{readout_list}")
    lines.append(f"The uploaded pooling strategy comprizes {n_compounds} samples in {n_pools} pools.")

    if n_compounds_f < 2:
        if n_compounds_f == 1:
            lines.append("A single positive sample was found:")
            lines.append(f"Sample: {decoded[0]}")
        else:
            lines.append("We found no matches for the given parameters, check your input or try increasing the differentiate value.")
        return "\n".join(lines)

    if decoded_type == "putative" and decoded_set is not None and len(decoded_set) > 0:
        if putative_reason == "max_combs":
            lines.append("Putative positive samples were identified, but the app does not have the computational power to attempt to decode the exact combination.")
            lines.append("Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.")
        else:
            lines.append("Putative positive samples were identified, but the exact combination could not be pinpointed.")
            lines.append("Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.")
        lines.append(f"There are up to {min([diff_deco, len(decoded_set)])} positive samples among the following samples: {decoded_set}.")
        return "\n".join(lines)

    if decoded_type == "error":
        lines.append("We found no matches for the given parameters, check your input or try increasing the differentiate value.")
    elif decoded_type == "unique":
        if isinstance(decoded, list) and len(decoded) == 1 and isinstance(decoded[0], list) and len(decoded[0]) == 1:
            lines.append("A single positive sample was found:")
            lines.append(f"Sample: {decoded[0][0]}")
        else:
            lines.append("A single possible combination of positive samples was found. The positive samples are:")
            combo = decoded[0] if isinstance(decoded, list) and len(decoded) == 1 else decoded
            lines.append(f"Samples: {', '.join(map(str, combo))}")
    else:
        lines.append(f"{len(decoded)} possible combinations of positive samples were found. The possible combinations are:")
        for i, deco in enumerate(decoded):
            if i != 0:
                lines.append("or")
            lines.append(f"Samples: {deco}")

    return "\n".join(lines)



def decode_continuous_lasso(readout_arr, WA, diff=None, min_signal=None, diluting=True, alpha=0.05, l1_ratio=0.5):
    try:
        n_compounds, n_pools = WA.shape
        readout_arr = np.array(readout_arr, dtype=float).flatten()
        
        if len(readout_arr) != n_pools:
            return "error", f"Readout length {len(readout_arr)} doesn't match pool count {n_pools}", None, None, n_compounds, None
        
        # Feature matrix: X = WA.T (n_pools, n_compounds)
        X = WA.T.astype(float)
        y = readout_arr
        
        # Use differentiate (diff) to control sparsity: lower alpha for more sparsity
        # diff represents max number of active compounds, so use it to scale alpha
        if diff is not None and diff > 0:
            # Scale alpha based on sparsity constraint
            # More compounds active -> higher alpha (less sparsity)
            alpha_scaled = alpha * (diff / n_compounds)
        else:
            alpha_scaled = alpha
        
        model = Lasso(alpha=alpha_scaled, positive=True, max_iter=10000)
        model.fit(X, y)
        # Get the optimized coefficients
        
        return "continuous",  model.coef_
    
    except Exception as e:
        return "error", f"Error in decode_continuous: {str(e)}"


def decode_continuous_elastic(readout_arr, WA, diff=None, min_signal=None, diluting=True, alpha=1.0, l1_ratio=0.5):
    try:
        n_compounds, n_pools = WA.shape
        readout_arr = np.array(readout_arr, dtype=float).flatten()
        
        if len(readout_arr) != n_pools:
            return "error", f"Readout length {len(readout_arr)} doesn't match pool count {n_pools}", None, None, n_compounds, None
        
        # Feature matrix: X = WA.T (n_pools, n_compounds)
        X = WA.T.astype(float)
        y = readout_arr
        model = ElasticNet(alpha=alpha,l1_ratio=l1_ratio, positive=True, max_iter=10000)
        model.fit(X, y)
        # Get the optimized coefficients
        
        return "continuous",  model.coef_
    
    except Exception as e:
        return "error", f"Error in decode_continuous: {str(e)}"

