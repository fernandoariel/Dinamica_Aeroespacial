import numpy as np 

def zero_to_360(x):
    if x>= 360:
        x=x-np.fix(x/360)*360
    elif x<0:
        x=x-(np.fix(x/360)-1)*360

    y=x
    return y      
