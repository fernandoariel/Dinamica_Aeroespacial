import numpy as np 

def J0(year, month, day):
    J0=367*year-np.fix(7*(year+np.fix((month+9)/12))/4)+np.fix(275*month/9)+day+1721013.5

    return J0

