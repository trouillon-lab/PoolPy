import numpy as np
import matplotlib.pyplot as plt
import itertools



#Coumpound counter start from 1
# Helper function for binary translation
def IntegerToBinaryTF(num: int, ls_bn: list)-> list:
    if num >= 2:
        ls_bn=IntegerToBinaryTF(num // 2, ls_bn)
    ls_bn.append(num % 2==1)
    return(ls_bn)


# Function to assigne each compound to its wells
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


# Function to be called by user to create the compound-wells assignment matrix
# This functions also identifies the minimum nuber of wells needed for the compunds and level of detail (differentiate) selected
def assing_wells_L(n_compounds:int, differentiate=1) -> np.array:

    if differentiate==1:
        n_wells=int(np.ceil(np.log2(n_compounds +1)))
    if differentiate==2:
        tentative=int(np.floor(np.sqrt(6*n_compounds))) #empirical evidence suggest almost this scaling, mathematical proof mught arrive later
        for NW in [tentative-1,tentative,tentative+1]:
            if get_ncomp_from_nwells(NW, differentiatie=2)>=n_compounds:
                n_wells=NW
                break
        
    well_assigner=np.zeros((n_compounds, n_wells))==1
    for i in range(n_compounds):
        well_assigner[i,:]=well_selecter(i+1, n_wells, differentiate)
    return(well_assigner)

def assign_wells_mat(n_compounds:int)->np.array:
    L1=np.ceil(np.sqrt(n_compounds))
    L2=L1-1 if L1*(L1-1)>=n_compounds else L1
    well_assigner=np.zeros((n_compounds, int(L1+L2)))==1
    for i in range(n_compounds):
        cp_id=[int(i//L1), int(L1+i % L1)]
        well_assigner[i,cp_id]=True
    return(well_assigner)


# Helper function to actually implement the experiment
def from_well_get_compuonds(well:int, well_assigner: np.array)->np.array :
    return(np.array(np.where(well_assigner[:,well-1]))[0]+1)

def from_compound_get_wells(compound: int, well_assigner: np.array)-> np.array:
    return(np.array(np.where(well_assigner[compound-1,:]))[1]+1)

# Consistency check

def is_consistent(well_assigner:np.array, differentiate:int) -> list:
    n_comp=well_assigner.shape[0]
    if differentiate==1:
        full_well_assigner=well_assigner.copy()
    if differentiate==2:
        n_perms=int(n_comp*(n_comp-1)/2+n_comp)
        full_well_assigner=np.zeros((n_perms, well_assigner.shape[1]))==1
        full_well_assigner[-n_comp:,:]=well_assigner.copy()
        for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),2)):
            k,l=j
            full_well_assigner[i,:]=well_assigner[k,:]+well_assigner[l,:]
    if differentiate==3:
        n_perms=int(n_comp*(n_comp-1)*(n_comp-2)/6+n_comp*(n_comp-1)/2+n_comp)
        n3=int(n_comp*(n_comp-1)*(n_comp-2)/6)
        full_well_assigner=np.zeros((n_perms, well_assigner.shape[1]))==1
        full_well_assigner[-n_comp:,:]=well_assigner.copy()
        for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),3)):
            k,l,m=j
            full_well_assigner[i,:]=well_assigner[k,:]+well_assigner[l,:]+well_assigner[m,:]
        for i, j in enumerate(itertools.combinations(np.arange(well_assigner.shape[0]),2)):
            k,l=j
            full_well_assigner[i+n3,:]=well_assigner[k,:]+well_assigner[l,:]
    if np.unique(full_well_assigner, axis=0).shape[0]<full_well_assigner.shape[0]:
        _, counts=np.unique(full_well_assigner, axis=0, return_counts=True)
        return(False, full_well_assigner, counts)
    elif np.unique(full_well_assigner, axis=0).shape[0]==full_well_assigner.shape[0]:
        return(True,full_well_assigner, counts)
    else:
        return("Something is fishy")
        

 

# Select number of compuonds and make the well assigner
NC=120
WA=assing_wells_L(NC,differentiate=2)

#check shape for consitency
WA.shape

#Print where a compound goes 
cp=np.random.randint(1,WA.shape[0]+1,1)
well_cp=from_compound_get_wells(cp,WA)
print('Compound', cp, 'needs to be in wells', well_cp)

# Get all compounds for a well
wl=np.random.randint(1,WA.shape[1]+1,1)
cp_well=from_well_get_compuonds(wl,WA)
print('In well', wl, 'need to go compunds', cp_well)

# Check answer ambiguity
tf, WAS, countss= is_consistent(WA, 2)

#Plot ambiguity
plt.hist(countss)
plt.xlabel('Possible ways to create an outcome readout')
plt.show()

pr1=np.sum(countss[countss>1])/np.sum(countss)
pr2=1-np.sum(countss[countss==1])/len(countss)
print('with this assignment method from',np.round(pr2*100) ,'to', np.round(pr1*100),'% of outcomes will lead to an ambiguous interpretation')

WAM=assign_wells_mat(120)

WAM.shape

tf, WAT, countsm= is_consistent(WAM, 2)

plt.hist(countsm)
plt.xlabel('Possible ways to create an outcome readout')
plt.show()

pr1=np.sum(countsm[countsm>1])/np.sum(countsm)
pr2=1-np.sum(countsm[countsm==1])/len(countsm)
print('with this assignment method from',np.round(pr2*100) ,'to', np.round(pr1*100),'% of outcomes will lead to an ambiguous interpretation')

cp=np.random.randint(1,WA.shape[0]+1,1)
well_cp=from_compound_get_wells(cp,WA)
print('Compound', cp, 'needs to be in wells', well_cp)
    

wl=np.random.randint(1,WA.shape[1]+1,1)
cp_well=from_well_get_compuonds(wl,WA)
print('In well', wl, 'need to go compunds', cp_well)
