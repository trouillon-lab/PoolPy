import numpy as np








def add_1(combinantions_dictionary, ND=5):
    N=combinantions_dictionary[1][-1]+1
    new_cd={1:np.append(combinantions_dictionary[1],N)}
    diff=2
    while diff<ND:
        new_part=np.vstack([combinantions_dictionary[diff-1].T,np.array([N]*len(combinantions_dictionary[diff-1]))])
        new_in=np.hstack([combinantions_dictionary[diff].T,new_part]).T
        new_cd.update({(diff):new_in})

        diff+=1






