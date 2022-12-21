import math as m
import numpy as np



def rotate (r,Vek):
    sp = m.sin((-r/360)*2*m.pi)
    cp = m.cos(-(r/360)*2*m.pi)
    r_curr = np.array([[1,0,0],[0,cp,sp],[0,-sp,cp]]) 
    RVek=r_curr.dot(Vek)
   
    return RVek
                
    