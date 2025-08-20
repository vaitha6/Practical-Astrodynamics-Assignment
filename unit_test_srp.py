'''
The concept used was to try to first calculate the angle between the vector from the Sun to the satellite and the Earth to the satellite, and use that angle
to find the vector from the Earth to Point P using trignometric relations. This concept was attempted using the given code but a math domain error kept
occuring when using the math.acos() function. This method was then scrapped and the formulae provided in the lecture slides were used.
'''


from math import *
from sgp4.api import *
import numpy as np
import matplotlib.pyplot as plt
from sgp4.conveniences import sat_epoch_datetime
from datetime import datetime,timedelta
import time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_sun, get_moon


mu=398600.4415E+9
mu_sun=1.32712440042E+20
mu_moon=4.9028695E+12
R_e=6378.1363E+3
R_sun=6.96E+8
G=6.667E-11
M=6E+24
J2=1.083E-3
J22=1.8155628E-6
lambda_22=-14.9287*(pi/180)
rho_0=1.2
H=8E+3
B=0.01
omega_e=7292115e-11
r_sun=np.array([-2.141664114391220E+10, -1.336901089766156E+11, -5.795020557937574E+10])
r_moon=np.array([-3.316594688659666E+08, 1.896187710272793E+08, 1.235349432352076E+08])
B_r=0.02
W=1361
c=299792458

out3=np.loadtxt('out3.dat',int)

r_es=np.array([out3[:,1],out3[:,2],out3[:,3]])

for i in range(7000):
    r_ss=r_es[:,i]-r_sun
    
    r_es_mag=sqrt((r_es[0,i]**2)+(r_es[1,i]**2)+(r_es[2,i]**2))
    r_ss_mag=sqrt((r_ss[0]**2)+(r_ss[1]**2)+(r_ss[2]**2))
    
    r_es_u=r_es/r_es_mag
    r_ss_u=r_ss/r_ss_mag
    
    r_es_u_mag=sqrt((r_es_u[0]**2)+(r_es_u[1]**2)+(r_es_u[2]**2))
    r_ss_u_mag=sqrt((r_ss_u[0]**2)+(r_ss_u[1]**2)+(r_ss_u[2]**2))
    theta=acos(r_ss_u_mag/r_es_u_mag)
    
    r_ep=r_es[:,i]*sin(theta)
    r_ep_mag=sqrt((r_ep[0]**2)+(r_ep[1]**2)+(r_ep[2]**2))
    h_g=r_ep_mag-R_e
    
    r_ps=r_es[:,i]-r_ep
    r_ps_mag=sqrt((r_ps[0]**2)+(r_ps[1]**2)+(r_ps[2]**2))
    R_p=(r_ps_mag*R_sun)/r_ss_mag
    
    shadow_func=h_g/R_p
