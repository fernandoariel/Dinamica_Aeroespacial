import numpy as np 
#Programa para encontrar los parametros orbitales partiendo de
# un vector de estado R y V conocidos
#################################
def orbital_elements(R,V,mu):
    eps=1E-10    
    r= np.linalg.norm(R)
    v= np.linalg.norm(V)
    v_r= np.dot(V, R)/r
    h_vector=np.cross(R, V)
    h= np.linalg.norm(h_vector)


    inc=np.arccos(h_vector[2]/h)

    N=np.cross(np.array([0, 0, 1]), h_vector)
    n=np.linalg.norm(N)
    if n != 0:
        RA=np.arccos(N[0]/n)
        if N[1]<0:
            RA=2*np.pi-RA
    else:
        RA=0

    E=(1/mu)*((v**2-mu/r)*R - r*v_r*V)
    e=np.linalg.norm(E)

    if n != 0:
        if e>eps:
            w=np.arccos(np.dot(N, E)/(n*e))
            if E[2]<0:
                w=2*np.pi - w
        else:
            w = 0
    else:
        w=0


    if e > eps:
        TA= np.arccos(np.dot(E, R)/(e*r))    
        if v_r < 0:
            TA= 2*np.pi-TA      
    else:
        cp=np.cross(N,R)
        if cp[2] >= 0:
            TA=np.arccos(np.dot(N,R)/(n*r))
        else:
            TA=2*np.pi - np.arccos(np.dot(N, R)/(n*r))

    a=h**2/(mu*(1-e**2))
    deg=180/np.pi

    return  (h, e, RA*deg, inc*deg, w*deg, TA*deg, a )

######################
# input

mu=398600
R=np.array([-6045, -3490, 2500])   
V=np.array([-3.457, 6.618, 2.533 ])

#output

h, e, RAA, i, w, theta, a = (orbital_elements(R, V, mu))

print('El momento angular es: %.4f km2/seg '%h)
print('La excentricidad de la orbita es: %.4f'%e)
print('El angulo de ascenscion recta del Nodo Ascendentees: %.4f deg'%RAA)
print('La inclinacion de la orbita es: %.4f deg'%i)
print('El argumento del perigeo es: %.4f deg'%w)
print('La anomalia verdadera es: %.4f deg'%theta)
print('El semieje mayor es: %.4f km'%a)