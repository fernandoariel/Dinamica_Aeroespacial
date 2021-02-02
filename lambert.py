import numpy as np
from coe_from_sv import coe_from_sv
def stumpS(z):
    if z>0:
        s=(np.sqrt(z)-np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    elif z<0:
        s=(np.sinh(np.sqrt(-z))-np.sqrt(-z))/(np.sqrt(-z))**3
    else:
        s=1/6
    return s    


def stumpC(z):
    if z>0:
        c=(1-np.cos(np.sqrt(z)))/z
    elif z<0:
        c=(np.cosh(np.sqrt(-z))-1)/(-z)
    else:
        c=1/2
    return c    


def lambert(R1, R2, t):
    mu=1.327124E11
    r1=np.linalg.norm(R1)
    r2=np.linalg.norm(R2)
    cl2=np.cross(R1, R2)
    theta= np.arccos(np.dot(R1,R2)/r1/r2)

    if cl2[2] <= 0:
        theta=2*np.pi - theta
    A= np.sin(theta)*np.sqrt((r1*r2)/(1-np.cos(theta)))

    def C(z):
        dum= stumpC(z)
        return dum
    def S(z):
        dum= stumpS(z)
        return dum  
    def y(z):
        dum=r1 + r2 + A*(z*S(z)-1)/np.sqrt(C(z))
        return dum
   
    def F(z,t):
        dum= (y(z)/C(z))**1.5*S(z) + A*np.sqrt(y(z))-np.sqrt(mu)*t  
        return dum
            
    def dFdz(z):
        if z == 0:
            dum= np.sqrt(2)/40*y(0)**1.5 + A/8*(np.sqrt(y(0))+A*np.sqrt(1/2/y(0)))
        else:
            dum= (y(z)/C(z))**1.5*(1/2/z*(C(z)-3*S(z)/2/C(z))+3*S(z)**2/4/C(z))+A/8*(3*S(z)/C(z)*np.sqrt(y(z)) + A*np.sqrt(C(z)/y(z)))
        return dum
    z=0
    while F(z, t)<0:
        z=z+0.1

    tol=1E-8
    nmax= 5000

    ratio=1
    n=0
    while np.abs(ratio)>tol and n<=nmax:
        n=n+1
        ratio= F(z,t)/dFdz(z)
        z=z-ratio    

    if n>= nmax:
        print('number of itirations exceeds')

    f=1-y(z)/r1

    g=A*np.sqrt(y(z)/mu) 
    gdot= 1 - y(z)/r2
    V1= 1/g*(R2-f*R1)
    V2= 1/g*(gdot*R2-R1)

    return V1,V2

