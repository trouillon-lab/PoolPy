from Functions import *



####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
#################################################################################################### DECODE ########################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

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
            idxs = np.all(outcome == full_well_assigner, axis=1)
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
                    idxs=np.prod(outcome==full_well_assigner, axis=1)
                    outcome_dict.update({tuple_to_str(tuple(outcome)):list(itertools.compress(sc_list,idxs))})
                full_od.update({diff:outcome_dict})

            else:
                idxs=np.prod(readout==full_well_assigner, axis=1)
            full_od.update({diff:list(itertools.compress(sc_list,idxs))})
        return full_od
    
def decode_sweep(dir_scramblers, dir_WAs, differentiate:int,
            max_differentiate,
            start=50, stop=150, step=10,
             **kwargs) -> list:
    N=start
    while N<stop:
        Npath=os.path.join(dir_WAs,'N_'+str(N))
        diff=1
        if max_differentiate>=1:

            diff=1
            if 'differentiate' in kwargs.keys():
                del kwargs['differentiate']
            while diff<=max_differentiate:
                start_time = time.time()
                dpath=os.path.join(Npath,'diff_'+str(diff))
                if diff==1:
                    scrambler={1:np.arange(N)}
                    WApath=os.path.join(dpath,'WAs')
                    filenames = next(os.walk(WApath), (None, None, []))[2]
                    for fname in filenames:
                        fdir=os.path.join(WApath,fname)
                        WA=np.genfromtxt(fdir, delimiter=",")
                        method=re.sub('^WA_', '', fname)
                        method=re.sub('_.*$', '', method)
                        dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 
                        scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                        decpath=os.path.join(dpath,'decoders')
                        if not os.path.exists(decpath):
                            os.makedirs(decpath)
                        decname=os.path.join(decpath, 'decoder_'+method+'.json')
                        #print(dict_decode)
                        json.dump( dict_decode, open(decname, 'w' ) )

                else:
                    this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                    if os.path.isfile(this_sc_file):
                        this_scrambler=np.load(this_sc_file)['sc']
                        scrambler.update({diff:this_scrambler})
                        dpath=os.path.join(Npath,'diff_'+str(diff))
                        WApath=os.path.join(dpath,'WAs')
                        filenames = next(os.walk(WApath), (None, None, []))[2]
                        for fname in filenames:
                            fdir=os.path.join(WApath,fname)
                            WA=np.genfromtxt(fdir, delimiter=",")
                            method=re.sub('^WA_', '', fname)
                            method=re.sub('_.*$', '', method)
                            dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 

                            scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                            decpath=os.path.join(dpath,'decoders')
                            if not os.path.exists(decpath):
                                os.makedirs(decpath)
                            decname=os.path.join(decpath, 'decoder_'+method+'.json')

                            json.dump( dict_decode, open(decname, 'w' ) )


                DTS=np.round((time.time() - start_time),2)
                DTD=DTS//86400
                DTH=DTS//3600-DTD*24
                DTM=DTS//60-DTH*60-DTD*24*60
                DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
                
                
                print('\n')
                print('----------------------------------------------------------------------------------------------------------') 
                print("%s days %s hours %s minutes and %s seconds required for N= %s and differentiate %s" % 
                        (DTD, DTH, DTM, DTS, N, diff))
                print('----------------------------------------------------------------------------------------------------------') 



                diff+=1

        else:
            diff=1    
            while diff<=differentiate:  
                if diff==1:
                        scrambler={1:np.arange(N)}
                        WApath=os.path.join(dpath,'WAs')
                        filenames = next(os.walk(WApath), (None, None, []))[2]
                        for fname in filenames:
                            fdir=os.path.join(WApath,fname)
                            WA=np.genfromtxt(fdir, delimiter=",")
                            method=re.sub('^WA_', '', fname)
                            method=re.sub('_.*$', '', method)
                            dict_decode=decode_precomp(well_assigner=WA, differentiate=diff, 
                            scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                            decname=os.path.join(dpath, 'decoder_'+method+'.json')
                            json.dump( dict_decode, open(decname, 'w' ) )

                else:
                    this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                    this_scrambler=np.load(this_sc_file)['sc']
                    scrambler.update({diff:this_scrambler})
                    dpath=os.path.join(Npath,'diff_'+str(diff))
                    WApath=os.path.join(dpath,'WAs')
                    filenames = next(os.walk(WApath), (None, None, []))[2]
                    for fname in filenames:
                        fdir=os.path.join(WApath,fname)
                        WA=np.genfromtxt(fdir, delimiter=",")
                        decode_precomp(well_assigner=WA, differentiate=diff, 
                        scrambler=scrambler, readout=np.nan, max_differentiate=-1, sweep=True, **kwargs)
                        decname=os.path.join(dpath, 'decoder_'+method+'.json')
                        json.dump( dict_decode, open(decname, 'w' ) )
 


        N+=step




def decode_single( WA, readout:np.ndarray, 
                  differentiate, scrambler=True,
            dir_scramblers=False, **kwargs) -> list:
        print('Function not yet ready, come back later')
        return 'Function not yet ready, come back later'
        if scrambler==True:
            N=WA.shape[0]
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        elif not scrambler==False:
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        
        elif not dir_scramblers==False:
            N=WA.shape[0]
            if differentiate==1:
                scrambler={1:np.arange(N)}
            return decode_precomp(well_assigner=WA, readout=readout, 
                                  differentiate=differentiate, scrambler=scrambler, max_differentiate=-1)
        

def str_to_tuple(string, delimiter='-'):
    return tuple(string.split(delimiter))

def tuple_to_str(tuple_type, delimiter='-'):
    return delimiter.join(map(str,tuple_type))
 
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
##################################################################################################### RANDOM ######################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################    



def get_max_C(n_compounds, max_compounds):
    return int(n_compounds/2) if max_compounds==0 else max_compounds
def get_min_C(n_compounds, MC):
    return int(np.sqrt(n_compounds)) if int(np.sqrt(n_compounds))<MC else int(MC/2)

def get_max_W(n_compounds):
    return int(np.log2(n_compounds))
def get_min_W(n_compounds):
    return int(2*np.sqrt(n_compounds))


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
            full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
    _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
    if len(counts)<full_well_assigner.shape[0]:
        return(False, full_well_assigner, counts)
    elif len(counts)==full_well_assigner.shape[0]:
        return(True,full_well_assigner, counts)
    else:
        print("Something is fishy")
        return(-1)


def mean_metrics_precomp(well_assigner, differentiate, scrambler, **kwargs):
    BT=well_assigner.shape[1]
    _,_, counts= is_consistent_precomp(well_assigner, differentiate, scrambler) 
    ET=extra_tests(counts)   
    rounds=np.sum(counts>1)/np.sum(counts>0)+1
    p_check=np.round(np.sum(counts[counts>1])/np.sum(counts)*100)
    return BT+ET, ET,  rounds, p_check



def find_rand_params_precomp(n_compounds:int, n_compounds_per_well=0, n_wells=0, guesses=0, 
                     max_compounds=0, max_redundancy=4, min_redundancy=1,**kwargs):
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

    arr_wells=np.arange(mw,MW)
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
                
            if comp*wells>max_redundancy*n_compounds*np.log2(n_compounds) or comp*wells<min_redundancy*n_compounds*np.log2(n_compounds): continue 
            WA_tmp, mean_exp, p_check=evaluate_rand_design(n_compounds=n_compounds,n_wells=wells,n_compounds_per_well=comp, return_me=True, guesses=guesses, **kwargs)
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


def evaluate_rand_design(n_compounds:int,  differentiate:int,scrambler:dict, n_compounds_per_well=0, 
                        n_wells=0, guesses=0,  return_me=False, **kwargs):
    min_tests=np.inf
    second_axis=np.tile(np.arange(n_wells),n_compounds_per_well).reshape(n_compounds_per_well,-1)
    for i in range(guesses):
        #print(i)
        idt=np.random.randint(0,n_compounds,size=(n_compounds_per_well,n_wells) )
        well_assigner=np.zeros((n_compounds,n_wells))==1
        well_assigner[idt, second_axis]=True
        if guesses==1 and False:
            if return_me:
                mean_exp, _, _, p_check= mean_metrics_precomp(well_assigner=well_assigner, 
                                                            differentiate=differentiate,scrambler=scrambler,**kwargs)
                return well_assigner, mean_exp, p_check
            return well_assigner
        mean_exp, _, _, p_check= mean_metrics_precomp(well_assigner=well_assigner,
                                                    differentiate=differentiate, scrambler=scrambler, **kwargs)
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


def assign_wells_random_precomp(n_compounds:int,  differentiate:int,scrambler:dict, n_compounds_per_well=0, 
                        n_wells=0, guesses=0, Evaluate=False, return_me=False, **kwargs)->np.array:

    _,_, min_tests, WA_rand, p_check=find_rand_params_precomp(n_compounds=n_compounds, differentiate=differentiate, 
                                 n_compounds_per_well=n_compounds_per_well, n_wells=n_wells, 
                                 guesses=guesses,scrambler=scrambler, **kwargs)
    if return_me:
        return WA_rand,  min_tests, p_check
    
    return WA_rand

def check_Rand_in_WApath(WApath):
    try:
        # List all files in the directory
        files = os.listdir(WApath)
        # Check if any file starts with 'WA_Random_N_'
        for file in files:
            if file.startswith('WA_Random_N_'):
                return True
        return False
    except Exception as e:
        return str(e)



def rand_sweep_diff(n_compounds, max_diff, dir_scramblers, Npath, **kwargs):
    if 'differentiate' in kwargs.keys():
        del kwargs['differentiate']
    N=n_compounds
    
    if max_diff>1:
        
        dpath=os.path.join(Npath,'diff_'+str(max_diff))
        WApath=os.path.join(dpath,'WAs')
        if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
            return
        
        
        for di in range(max_diff):
            start_time = time.time()
            diff=di+1
            if diff==1:
                dpath=os.path.join(Npath,'diff_'+str(diff))
                WApath=os.path.join(dpath,'WAs')
                scrambler={1:np.arange(n_compounds)}

                if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
                    continue

                
                WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                                differentiate=diff,scrambler=scrambler, return_me=True, **kwargs )
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

                #.append(['Random', min_tests, np.max(np.sum(WA_rand, axis=0)), WA_rand.shape[0], int(perc_check),  extra_exp,1+perc_check/100])
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
                np.savetxt(thisfile, WA_rand.astype(bool), delimiter=",")


            else:
                dpath=os.path.join(Npath,'diff_'+str(diff))
                WApath=os.path.join(dpath,'WAs')
                this_sc_file=os.path.join(dir_scramblers, 'N_'+str(N),  'N_'+str(N)+'_diff_'+str(diff)+'.npz')
                this_scrambler=np.load(this_sc_file)['sc']
                scrambler.update({diff:this_scrambler})

                if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
                    continue


                WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                                differentiate=diff,scrambler=scrambler, return_me=True, **kwargs  )
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


                #.append(['Random', min_tests, np.max(np.sum(WA_rand, axis=0)), WA_rand.shape[0], int(perc_check),  extra_exp,1+perc_check/100])
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
                np.savetxt(thisfile, WA_rand.astype(bool), delimiter=",")

            DTS=np.round((time.time() - start_time),2)
            DTD=DTS//86400
            DTH=DTS//3600-DTD*24
            DTM=DTS//60-DTH*60-DTD*24*60
            DTS=np.round(DTS-(DTM+DTH*60+DTD*24*60)*60,2)
            print("%s days %s hours %s minutes and %s seconds required for N= %s and differentiate %s" % 
                  (DTD, DTH, DTM, DTS, n_compounds, diff))
            print('----------------------------------------------------------------------------------------------------------') 

    elif max_diff==1:

        start_time = time.time()
        diff=1
        dpath=os.path.join(Npath,'diff_'+str(diff))
        WApath=os.path.join(dpath,'WAs')

        if not kwargs['overwrite'] and check_Rand_in_WApath(WApath):
            return

        scrambler={1:np.arange(n_compounds)}
        WA_rand,  min_tests, perc_check=assign_wells_random_precomp(n_compounds=n_compounds, 
                                                        differentiate=diff,scrambler=scrambler, return_me=True, **kwargs )
        extra_exp=WA_rand.shape[1]+min_tests
        #.append(['Random', min_tests, np.max(np.sum(WA_rand, axis=0)), WA_rand.shape[0], int(perc_check),  extra_exp,1+perc_check/100])
        
        
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
        np.savetxt(thisfile, WA_rand.astype(bool), delimiter=",")


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



####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
################################################################################################## DETERMINISTIC ###################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################
####################################################################################################################################################################################################################













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
    timeit = kwargs.get('timeit', False)
    
    # Generate multidim methods dynamically
    multidim_methods = []
    for i in np.arange(2, int(np.ceil(np.log(n_compounds)/np.log(2)))):
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
    for diff in range(1, max_diff + 1):
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
                    np.savetxt(dst_file, WA_list[idx].astype(bool), delimiter=",")
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
        for method in ['STD', 'Chinese trick']:
            dst_file = get_wa_filename(save_dir, n_compounds, diff, method)
            
            if not os.path.exists(dst_file):
                # Compute only if file doesn't exist
                if method == 'STD':
                    WA = assign_wells_STD(**current_kwargs)
                else:  # Chinese trick
                    WA = assign_wells_chinese(**current_kwargs)
                
                np.savetxt(dst_file, WA.astype(bool), delimiter=",")
                if timeit:
                    print(f"Computed {method} for diff={diff}")
            elif timeit:
                print(f"Skipping {method} for diff={diff} (already exists)")

def make_all_deterministic_WAs(start=50, stop=150, step=10, **kwargs):
    """
    Main loop to process all n_compounds values
    """
    current = start
    while current < stop:
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

# Argument parsing and main execution remains the same as in your original code
# ... [rest of your argument parsing and main call] ...




def make_all_deterministic_WAs_old(start=50, stop=150, step=10, **kwargs):
    
    current=start

    while current<stop:
        time0=time.time()
        if kwargs['timeit']:
            print(current)
        full_deterministic_WAS(n_compounds=current, **kwargs)
        current=current+step
        if kwargs['timeit']:
            print("segment time: %s seconds" % np.round(time.time() - time0, 1))


def full_deterministic_WAS(**kwargs):
    #methods=['matrix', 'random', 'STD', 'Chinese trick']
    # matrix assignment
    
    kwargs['return_wa']=True

    WA_list=[]
    multi=[]
    for i in np.arange(2,int(np.ceil(np.log(kwargs['n_compounds'])/np.log(2)))):
        if i>kwargs['max_dims']:
            continue
        WA_mul=assign_wells_multidim(n_dims=i, **kwargs)
        WA_list.append(WA_mul)
        multi.append('multidim-'+str(i))
        
    methods=multi.copy()

    WA_mat=assign_wells_mat(**kwargs)

    WA_list.append(WA_mat)
    methods.append('Matrix')

    WA_bin=assign_wells_bin(**kwargs)
    methods.append('Binary')
    WA_list.append(WA_bin)

    if kwargs['max_diff']>1:
        for diffo in np.arange(kwargs['max_diff']):


            kwargs.update({'differentiate':int(diffo+1)})

            WA_listo=copy.deepcopy(WA_list)
            methodso=copy.deepcopy(methods)

            # STD asignment 
            WA_std=assign_wells_STD(**kwargs)

            

            # chinese trick assignment
            WA_chin=assign_wells_chinese(**kwargs)


            WA_listo.extend([ WA_std, WA_chin])
            methodso.extend(['STD', 'Chinese trick'])



            this_dir=os.path.join(kwargs['save_dir'],'N_'+str(kwargs['n_compounds']), 'diff_'+str(kwargs['differentiate']), 'WAs')

            if not os.path.exists(this_dir):
                os.makedirs(this_dir)

            for method, WA in zip(methodso, WA_listo):
                thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_'+str(kwargs['differentiate'])+'.csv')
                np.savetxt(thisfile, WA.astype(bool), delimiter=",")
    
    else:
        


        WA_listo=copy.deepcopy(WA_list)
        methodso=copy.deepcopy(methods)

        # STD asignment 
        WA_std=assign_wells_STD(**kwargs)

        

        # chinese trick assignment
        WA_chin=assign_wells_chinese(**kwargs)


        WA_listo.extend([ WA_std, WA_chin])
        methodso.extend(['STD', 'Chinese trick'])



        this_dir=os.path.join(kwargs['save_dir'],'N_'+str(kwargs['n_compounds']), 'diff_1', 'WAs')

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        for method, WA in zip(methods, WA_list):
            thisfile=os.path.join(this_dir,'WA_'+ method+'_N_'+str(kwargs['n_compounds'])+'_diff_1.csv')
            np.savetxt(thisfile, WA.astype(bool), delimiter=",")


