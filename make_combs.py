import numpy as np







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






