import numpy as np 

######################
#Este codigo permite obtener el vector de estado R y V de un satelite que orbita la tierra luego
# de un determinado tiempo t. Es necesario el vector de estado R_0 y V_0 en el tiempo inicial
# El procedimiento se tomo del libro Curtis-Orbital Mechanics
######################
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

def f_and_g(x,t,ro,a):
    global mu
    z=a*x**2
    #ec. 3.69a
    f=1-(x**2/ro)*stumpC(z)
    #ec. 3.69b
    g=t-(1/np.sqrt(mu))*(x**3)*stumpS(z)

    return (f,g)

def fDot_and_gDot(x,r,ro,a):
    global mu
    z=a*x**2
    #ec. 3.69c
    fdot=(np.sqrt(mu)/(r*ro))*(z*stumpS(z)-1)*x
    #ec. 3.69d
    gdot=1-(x**2/r)*stumpC(z)
    return (fdot,gdot)


def kepler_u(dt,r0,vr0, a):
    global mu
    error=1E-8
    nmax=1000
    # calculo de la anomalia universal
    x=np.sqrt(mu)*np.abs(a)*dt

    n=0
    ratio=1
    while np.abs(ratio)>error or n<=nmax:
        n=n+1
        C=stumpC(a*x**2)
        S=stumpS(a*x**2)
        F=r0*vr0/np.sqrt(mu)*x**2*C + (1-a*r0)*x**3*S + r0*x -np.sqrt(mu)*dt
        dFdx=r0*vr0/np.sqrt(mu)*x*(1-a*x**2*S) + (1-a*r0)*x**2*C + r0
        ratio=F/dFdx
        x=x-ratio
        return x

def rv_from_r0_v0(R0, V0, t):
    global mu
    #magnitudes de R0 y V0
    r_0=np.linalg.norm(R0)  
    v_0=np.linalg.norm(V0)
    tf=t
    #velocidad radial inicial
    v_r_0=np.dot(R0,V0)/r_0
    #inversa del semiejemayor (partiendo de la ec. de la energia)
    alpha=2/r_0 - (v_0**2)/mu
    #calculo de la anomalia uniersal
    x=kepler_u(tf, r_0, v_r_0, alpha)
    #calculo de las funciones f y g
    f,g=f_and_g(x,tf,r_0,alpha)
    #calculo del vector posicion final
    R=f*R0+g*V0
    #calculo de la magnitud de R
    r=np.linalg.norm(R)
    #calculo de las derivadas de f y g
    fdot, gdot=fDot_and_gDot(x,r,r_0, alpha)
    #calculo del vector celocidad final
    V= fdot*R0+gdot*V0

    return (R,V)

#datos del problema

mu=398600 #parametro gravitacional
R_0=np.array([1600, 5310, 3800]) #vector R_0 inicial
V_0=np.array([-7.350, 0.4600, 2.470]) #vector V_0 inicial
t=3200 #tiempo final en seg
R,V=rv_from_r0_v0(R_0, V_0, t)

print('vector posición: ', R)
print('vector velocidad: ', V)
