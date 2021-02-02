import numpy as np 
from J0 import J0
from zero_to_360 import zero_to_360

def LST(y, m, d, ut, EL):
    j0= J0(y,m,d)
    j=(j0 - 2451545)/36525
    g0 = 100.4606184 + 36000.77004*j +0.000387933*j**2 - (2.583E-8)*j**3
    g0=zero_to_360(g0)
    gst = g0 + 360.98564724*(ut/24)
    
    lst_1= gst + EL
    lst = lst_1 - 360*np.fix(lst_1/360)

    return lst, j0



degrees= -64
minutes= -10
seconds= -51.78

#Date
year= 2020
month= 12
day= 11

#universal time

hour= 8
minute= 6
second= 6

#convert negative (west) longitude to east longitude

if degrees < 0:
    degrees= degrees + 360
#express the longitudes as decimal numbers
EL= degrees + minutes/60 + seconds/3600
WL= 360 - EL


#express universal time as a decimal number
ut=hour + minute/60 + second/3600    

lst, j0= LST(year, month, day, ut, EL)



print(lst)
print(j0)


