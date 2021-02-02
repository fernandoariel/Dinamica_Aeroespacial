import numpy as np 
from J0 import J0
from planetary_elements import planetary_elements
from zero_to_360 import zero_to_360
from kepler_E import kepler_E
from sv_from_coe import sv_from_coe

def planet_elements_and_sv(planet_id, year, month, day, hour, minute, second):
    mu=1.327124E11
    deg=np.pi/180

    #equation 5.47
    j0 = J0(year, month, day)

    ut = (hour + minute/60 +second/3600)/24

    #equation 5.47
    jd = j0 + ut

    J2000_coe, rates=planetary_elements(planet_id)

    t0=(jd-2451545)/36525

    elements=J2000_coe + rates*t0

    a= elements[0]
    e=elements[1]

    h=np.sqrt(mu*a*(1-e**2))

    incl=elements[2]
    RA=zero_to_360(elements[3])
    w_hat=zero_to_360(elements[4])
    L=zero_to_360(elements[5])
    w=zero_to_360(w_hat-RA)
    M=zero_to_360(L-w_hat)

    E=kepler_E(e, M*deg)
    
    TA=zero_to_360(2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))/deg)

    coe=np.array([h, e, RA, incl, w, TA, a, w_hat, L, M, E/deg])

    r, v =sv_from_coe(np.array([h, e, RA*deg, incl*deg, w*deg, TA*deg]), mu)

    return (coe, r, v, jd)

