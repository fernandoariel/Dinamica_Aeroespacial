import numpy as np 
#####################################
#  El siguiente código permite determinar las componentes del vector de estado en el 
#  sistema perifocal y en el sistema geocentrico ecuatorial partiendo de parametros orbitales
#  conocidos       
####################################
def vector_state(mu, e,i,omega,w,theta,h):
    V_1=np.array([np.cos(theta), np.sin(theta), 0])
    V_2=np.array([-np.sin(theta), e+np.cos(theta), 0])
    r_x=(h**2/mu)*(1/(1+e*np.cos(theta)))*V_1
    v_x=(mu/h)*V_2

    R3_raa= np.array([[np.cos(omega), np.sin(omega), 0], [-np.sin(omega), np.cos(omega), 0 ],[0, 0, 1]])
    R1_i= np.array([[1, 0, 0], [0, np.cos(i), np.sin(i)], [0, -np.sin(i), np.cos(i)]])
    R3_w= np.array([[np.cos(w), np.sin(w), 0], [-np.sin(w), np.cos(w), 0 ],[0, 0, 1]])

    Q=np.linalg.multi_dot([R3_w, R1_i, R3_raa])
    Q_px=np.transpose(Q)

    r=np.dot(Q_px, r_x)
    v=np.dot(Q_px, v_x)

    r=np.transpose(r)
    v=np.transpose(v)
    print('Las componentes del vector de estado en el sistema perifocal son: ')
    print('rx= ', r_x)
    print('vx= ', v_x)
    print('Las componentes del vector de estado en el sistema geocentrico ecuatorial son: ')
    print('r= ', r)
    print('v= ', v)
    




mu=398600 #constante km3/seg
e= 0.3 #excentricidad
A_p= 380 #altitud del perigeo (km)
i= 35 #angulo de inclinacion(º)
Omega=130 #RAAN (º)
w=115 # argumento del perigeo (°)
theta=0 #anomalia verdadera (°)
r_e=6378 # radio terrestre (km)
r_p=r_e+A_p #radio del perigeo
h=np.sqrt(r_p*mu*(1+e*np.cos(theta))) #momento angular
deg=np.pi/180  #conversion a radianes

vector_state(mu, e, i*deg, Omega*deg, w*deg, theta*deg, h)
