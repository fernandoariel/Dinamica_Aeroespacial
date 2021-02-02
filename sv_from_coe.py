import numpy as np 

def sv_from_coe(coe, mu):
    h=coe[0]
    e=coe[1]
    RA=coe[2]
    incl=coe[3]
    w=coe[4]
    TA=coe[5]

    rp=(h**2/mu)*(1/(1+e*np.cos(TA)))*np.array([np.cos(TA),np.sin(TA),0])

    vp=(mu/h)*np.array([-np.sin(TA), e + np.cos(TA), 0]) 

    R3_W=np.array([[np.cos(RA), np.sin(RA), 0],[-np.sin(RA), np.cos(RA), 0 ],[0, 0, 1]])

    R1_i=np.array([[1, 0, 0], [0, np.cos(incl), np.sin(incl)], [0, -np.sin(incl), np.cos(incl)]])

    R3_w=np.array([[np.cos(w), np.sin(w), 0],[-np.sin(w), np.cos(w), 0 ],[0, 0, 1]])

    Q_pX =np.linalg.multi_dot([R3_w,R1_i,R3_W])
    Q_pX=np.transpose(Q_pX)

    r=np.dot(Q_pX,rp)
    v=np.dot(Q_pX,vp)

    #r=np.transpose(r)
    #v=np.transpose(v)

    return (r, v)

