from Functions_wrapped import *
import matplotlib.pyplot as plt
from scipy import optimize as opt
from scipy import stats as sts


def ET_BG2_opt(N):
    return((1+(1-(1-p)**N)*N)/N)

def better_BG2_full(p):
    def ET_BG2_opt(N):
        return((1+(1-(1-p)**N)*N)/N)
    optimum=opt.minimize(ET_BG2_opt, np.sqrt(1/p))
    o=optimum['x'][0]
    o1=np.floor(o)
    o2=np.ceil(o)
    e1=ET_BG2_opt(o1)
    e2=ET_BG2_opt(o2)
    return [o1, e1] if e1<e2 else [o2, e2]

def better_BG2(p):
    def ET_BG2_opt(N):
        return((1+(1-(1-p)**N)*N)/N)
    optimum=opt.minimize(ET_BG2_opt, np.sqrt(1/p))
    o=optimum['x'][0]
    fv=np.round(o)
    e1=ET_BG2_opt(fv)
    return [fv, e1]

def better_BG2_D3(p):
    def ET_BG2_D3_opt(N):
        P0=1-(1-p)**(N[0]*N[1])
        P1=1-(1-p)**N[1]
        return((1+(P0*(N[0]+
                ((N[1]*N[0]*P1)/(1-(1-P1)**N[0]) ))) )/(N[0]*N[1]))
    base=np.power(1/p, 1/3)
    optimum=opt.minimize(ET_BG2_D3_opt, np.array([base,base]), bounds=[(1,None),(1,None)])
    o=np.array(optimum['x'])
    fv=np.round(o)
    e1=ET_BG2_D3_opt(fv)
    return [fv, e1]



def better_BG2_Gen(p, L):
    def ET_BG2_Gen_opt(N):
        Mm=np.cumprod(N)
        M=np.flip(np.cumprod(np.flip(N)))
        Mm=np.insert(Mm, 0, 1, axis=0)
        M=np.insert(M, len(M), 1, axis=0)
        N=np.insert(N, 0, 1, axis=0)
        psi=1-(1-p)**M
        return(((1+psi[0]*(np.sum(Mm[1:]*psi[:-1]/(1-(1-psi[:-1])**Mm[:-1]))) )/M[0]))
    base=np.power(1/p, 1/(L+1))
    optimum=opt.minimize(ET_BG2_Gen_opt, np.array([base]*L), bounds=[(1,None)]*L)
    o=np.array(optimum['x'])
    dl=0
    om=np.min(o)
    for i in o:
        if(i>10*om or i<1.5):
            dl+=1
    if dl>0:
        #print("Too many layers for optimal solution; optimum obtained with %d layers" %(L-dl))
        return(better_BG2_Gen(p, L-dl))
    fv=np.round(o)
    e1=ET_BG2_Gen_opt(fv)
    return [fv, e1]

def brute_better_BG2(p, ML):
    def ET_BG2_Gen_opt(N):
        Mm=np.cumprod(N)
        M=np.flip(np.cumprod(np.flip(N)))
        Mm=np.insert(Mm, 0, 1, axis=0)
        M=np.insert(M, len(M), 1, axis=0)
        N=np.insert(N, 0, 1, axis=0)
        psi=1-(1-p)**M
        return(((1+psi[0]*(np.sum(Mm[1:]*psi[:-1]/(1-(1-psi[:-1])**Mm[:-1]))) )/M[0]))
    e0=2
    for L in np.arange(1,ML+1):
        base=np.power(1/p, 1/(L+1))
        optimum=opt.minimize(ET_BG2_Gen_opt, np.array([base]*L), bounds=[(1,None)]*L)
        o=np.array(optimum['x'])
        fv=np.round(o)
        e1=ET_BG2_Gen_opt(fv)
        if e1<e0:
            e0=e1
            fv0=np.round(fv)
        else:
            return([fv0, e0])
    return([fv0, e0])
        



def better_BG2_Gen_MD(p, D):
    def ET_BG2_Gen_opt_MD(N):
        N=int(np.round(N))
        NS=N**D
        probs=[]
        rg=np.arange(N+1)
        for i in rg:
            probs.append(sts.binom.pmf(i,NS,p))
            if i==N:
                probs[-1]=sts.binom.sf(i,NS,p)
        ev=rg**D*probs
        ev[1]=0
        return((D*N+np.sum(ev))/(N**D))
    base=np.power(1/p, 1/(D+1))
    optimum=opt.minimize(ET_BG2_Gen_opt_MD, np.array(base), bounds=[(1,np.power(10/p, 1/D))])
    o=np.array(optimum['x'])
    fv=np.round(o)
    e1=ET_BG2_Gen_opt_MD(fv)
    return [fv, e1]

def ET_BG2_Gen_opt_brute_MD(N,D,p):
    N=int(np.round(N))
    NS=N**D
    probs=[]
    NV=[]
    NVV=[]
    rg=np.arange(N+1)
    for i in rg:
        probs.append(sts.binom.pmf(i,NS,p))
    probs=np.max([probs,[0]*len(probs)], axis=0)
    probs[-1]=1-np.sum(probs[:-1])
    probs=np.max([probs,[0]*len(probs)], axis=0)
    probs/=np.sum(probs)
    ev=rg**D*probs
    ev[1]=0
    return((D*N+np.sum(ev))/(N**D))

def brute_BG2_Gen_MD(p, D_max):
    ls_Ds=[]
    ls_lati=[]
    sulti=[]
    for D in np.arange(2,D_max+1):
        DI=False
        mino=2
        maxo=np.ceil(np.power(1/p, 1/D))+3
        lati=np.arange(mino,maxo)
        for lato in lati:
            sulto=ET_BG2_Gen_opt_brute_MD(lato, D,p)
            if sulto<2:
                sulti.append(sulto)
                ls_Ds.append(D)
                ls_lati.append(lato)
                DI=True
            else:
                break
        if not DI:
            loco=np.argmin(sulti)
            return [ls_Ds[loco], ls_lati[loco], sulti[loco]]
    loco=np.argmin(sulti)
    return [ls_Ds[loco], ls_lati[loco], sulti[loco]]