import numpy as np 
from planet_elements_and_sv import planet_elements_and_sv
from lambert import lambert
from coe_from_sv import coe_from_sv

def interplanetary(depart, arrive):
    mu=1.327124E11
    #depart
    planet_id_1= depart[0]
    year_1= depart[1]
    month_1= depart[2]
    day_1= depart[3]
    hour_1 = depart[4]
    minute_1 = depart[5]
    second_1= depart[6]
    #arrive
    planet_id_2= arrive[0]
    year_2= arrive[1]
    month_2= arrive[2]
    day_2= arrive[3]
    hour_2 = arrive[4]
    minute_2 = arrive[5]
    second_2= arrive[6]


    dum_1, Rp1, Vp1, jd1 = planet_elements_and_sv(planet_id_1, year_1, month_1, day_1, hour_1, minute_1, second_1)
    dum_2, Rp2, Vp2, jd2 = planet_elements_and_sv(planet_id_2, year_2, month_2, day_2, hour_2, minute_2, second_2)

    tof=(jd2-jd1)*24*3600

    #Pached conic asumption
    R1=Rp1
    R2=Rp2
    
    #velocity at departure and arrival assuming a prograde trajectory
    V1, V2 = lambert(R1, R2, tof)
    planet_1=np.array([Rp1, Vp1, jd1])
    planet_2=np.array([Rp2, Vp2, jd2])
    trajectory=np.array([V1, V2])

    return planet_1, planet_2, trajectory

mu=1.327124E11
deg=np.pi/180

#departure
planet_id_1= 2 #earth
year_1= 2005
month_1 = 12
day_1 = 1
hour_1= 0
minute_1 = 0
second_1 = 0

depart=np.array([planet_id_1, year_1, month_1, day_1, hour_1, minute_1, second_1])

#arrive
planet_id_2= 1 #venus
year_2= 2006
month_2 = 4
day_2 = 1
hour_2= 0
minute_2 = 0
second_2 = 0

arrive =np.array([planet_id_2, year_2, month_2, day_2, hour_2, minute_2, second_2])

planet_1, planet_2, trajectory = interplanetary(depart, arrive)

R_1= planet_1[0][0:3]
V_p_1= planet_1[1][0:3]
j_d_1= planet_1[2]

R_2= planet_2[0][0:3]
V_p_2= planet_2[1][0:3]
j_d_2= planet_2[2]

V_1 = trajectory[0]
V_2 = trajectory[1]

tof= j_d_2 - j_d_1

coe_1= coe_from_sv(R_1, V_1, mu)
coe_2= coe_from_sv(R_2, V_2, mu)

vinf_1 =V_1 - V_p_1
vinf_2 = V_2 - V_p_2

print('Departure from Earth \n ##################')
print('planet position vector: ', R_1)
print('magnitude: ', np.linalg.norm(R_1))
print('planet velocity: ', V_p_1)
print('magnitude: ', np.linalg.norm(V_p_1))
print('spacecraft velocity: ', V_1)
print('magnitude: ', np.linalg.norm(V_1))
print('v-infinity at departure: ', vinf_1)
print('magnitude: ', np.linalg.norm(vinf_1))
print('time of flight (days): ', tof)
print('###########################')
print('Arrival to Venus \n ###########################')
print('planet position vector: ', R_2)
print('magnitude: ', np.linalg.norm(R_2))
print('planet velocity: ', V_p_2)
print('magnitude: ', np.linalg.norm(V_p_2))
print('spacecraft velocity: ', V_2)
print('magnitude: ', np.linalg.norm(V_2))
print('v-infinity at arrival: ', vinf_2)
print('magnitude: ', np.linalg.norm(vinf_2))
print('time of flight (days): ', tof)

