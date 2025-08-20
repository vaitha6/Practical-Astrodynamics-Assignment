from math import *
from sgp4.api import *
import numpy as np
import matplotlib.pyplot as plt
from sgp4.conveniences import sat_epoch_datetime
from datetime import datetime,timedelta


mu=398600.4418E+9
R_e=6378.137E+3
G=6.667E-11
M=6E+24

sgp4_1='sgp4_size_5.dat'
sgp4_2='sgp4_size_1.dat'
integrator_1='euler_size_5.dat'
integrator_2='euler_size_1.dat'

sgp4_1_dat=np.loadtxt(sgp4_1, dtype=float)
sgp4_2_dat=np.loadtxt(sgp4_2, dtype=float)
integrator_1_dat=np.loadtxt(integrator_1, dtype=float)
integrator_2_dat=np.loadtxt(integrator_2, dtype=float)

dim_1=np.shape(sgp4_1_dat)
count_1=dim_1[0]

dim_2=np.shape(sgp4_2_dat)
count_2=dim_2[0]

diff_1=np.zeros([count_1,7])
diff_2=np.zeros(7)
local_err_1=np.zeros(count_1)
local_err_2=np.zeros(count_2)

time_1=np.zeros(count_1)
time_1=sgp4_1_dat[:,0]

for i in range(count_1):
    
    diff_1[i,:]=sgp4_1_dat[i,:]-integrator_1_dat[i,:]
    local_err_1[i]=sqrt((diff_1[i,1]**2)+(diff_1[i,2]**2)+(diff_1[i,3]**2)+(diff_1[i,4]**2)+(diff_1[i,5]**2)+(diff_1[i,6]**2))

time_2=np.zeros(count_2)
time_2=sgp4_2_dat[:,0]

for i in range(count_2):
    
    diff_2=sgp4_2_dat[i,:]-integrator_2_dat[i,:]
    local_err_2[i]=sqrt((diff_2[1]**2)+(diff_2[2]**2)+(diff_2[3]**2)+(diff_2[4]**2)+(diff_2[5]**2)+(diff_2[6]**2))

plt.plot(time_1,local_err_1)
plt.plot(time_2,local_err_2)
