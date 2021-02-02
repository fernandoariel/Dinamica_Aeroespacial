import numpy as np 
from coe_from_sv import coe_from_sv

#algoritmo para calcular el vector de estado a partir de los datos de la observacion
def rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H):
    global f, Re, wE, rad
    
    
    omega=  np.array([0, 0, wE])


    A=A*rad
    Adot= Adot
    a= a*rad
    adot=adot
    theta= theta*rad
    phi= phi*rad

       
    R=np.array([(Re/np.sqrt(1-(2*f-f**2)*(np.sin(phi))**2) + H)*np.cos(phi)*np.cos(theta),
               (Re/np.sqrt(1-(2*f-f**2)*(np.sin(phi))**2) + H)*np.cos(phi)*np.sin(theta), 
               (Re*(1-f)**2/np.sqrt(1-(2*f-f**2)*(np.sin(phi))**2)+H)*np.sin(phi)])

    Rdot=np.cross(omega, R)

    dec= np.arcsin(np.cos(A)*np.cos(phi)*np.cos(a) + np.sin(phi)*np.sin(a))

    h= np.arccos((np.cos(phi)*np.sin(a)-np.sin(phi)*np.cos(A)*np.cos(a))/np.cos(dec))

    if A>0 and A<np.pi:
        h= 2*np.pi - h
    
    RA= theta - h

    Rho=np.array([np.cos(dec)*np.cos(RA), np.cos(dec)*np.sin(RA), np.sin(dec)])

    r= R + rho*Rho

    decdot= (1/np.cos(dec))*(-Adot*np.cos(phi)*np.sin(A)*np.cos(a) + adot*(np.sin(phi)*np.cos(a) - np.cos(phi)*np.cos(A)*np.sin(a)))

    RAdot= wE + (Adot*np.cos(A)*np.cos(a) - adot*np.sin(A)*np.sin(a) + decdot*np.sin(A)*np.cos(a)*np.tan(dec))/(np.cos(phi)*np.sin(a) - np.sin(phi)*np.cos(A)*np.cos(a))

    Rhodot=np.array([-RAdot*np.sin(RA)*np.cos(dec) - decdot*np.cos(RA)*np.sin(dec), 
                    RAdot*np.cos(RA)*np.cos(dec) - decdot*np.sin(RA)*np.sin(dec),
                    decdot*np.cos(dec) ])

    v= Rdot + rhodot*Rho + rho*Rhodot

    return r, v, R, Rho, h, dec


rad=np.pi/180 #conversion a radianes
f=1/298.256421867 #factor de achatamiento de la tierra
Re=6378.13655 #Radio de la Tierra
wE= 7.292115E-5 # Velocidad angular de la tierra rad/s
mu=398600.4418 

#datos iniciales
rho=408
rhodot= 0
A= -150
Adot= 1E-3
a= 11
adot= 1E-3
theta= 137.8468
phi= -32.0967
H= 670 

r, v, R, Rho, h, dec =rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H )

coe= coe_from_sv(r, v, mu)



h=coe[0]
e=coe[1]
RA= coe[2]
incl= coe[3]
w= coe[4]
TA= coe[5]
a=coe[6]

rp=h**2/mu/(1+e)


print('Vector r: ', r)
print('Vector v: ',v)
print('Elementos orbitales:  ', coe)

