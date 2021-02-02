import numpy as np
def kepler_E(e, M):
    error=1E-8
    if M < np.pi:
        E= M + e/2
    else:
        E= M - e/2    

    ratio=1
    while np.abs(ratio)>error:
        ratio= (E -  e*np.sin(E)- M)/(1-e*np.cos(E))
        E=E-ratio
    return E 

