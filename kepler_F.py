import numpy as np 

def kepler_H(e, M):
    error=1E-8
    F=M
    ratio= 1
    while np.abs(ratio)  > error:
        ratio=(e*np.sinh(F)-F-M)/(e*np.cosh(F)-1)
        F=F-ratio
    return F

e=2.38
M=2.93

F=kepler_H(e, M)
print(F)
