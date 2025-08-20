import numpy as np
from math import *
import matplotlib.pyplot as plt


mu=398600.4418E+9
R_e=6378.137E+3
G=6.667E-11
M=6E+24

filename='out2.dat'
start=1000
stop=400000
step=5

count=int((stop-start)/step)

cart_data=np.loadtxt(filename, dtype=float)

x_init=cart_data[0,1]
y_init=cart_data[0,2]
z_init=cart_data[0,3]

x_dot_init=cart_data[0,4]
y_dot_init=cart_data[0,5]
z_dot_init=cart_data[0,6]

time=np.zeros(count)
time[0]=0

r=np.zeros([3,count])
r_mag=np.zeros(count)
r[:,0]=[x_init,y_init,z_init]
r_mag[0]=sqrt((x_init**2)+(y_init**2)+(z_init**2))

v=np.zeros([3,count])
v[:,0]=[x_dot_init,y_dot_init,z_dot_init]

acc=np.zeros([3,count])
acc[:,0]=-G*M*r[:,0]/(r_mag[0]**3)

y=np.zeros([6,count])
y_dot=np.zeros([6,count])

y[:,0]=np.hstack([r[:,0],v[:,0]])
y_dot[:,0]=np.hstack([v[:,0],acc[:,0]])


for i in range(count-1):
      
      y[:,i+1]=y[:,i]+(step*y_dot[:,i])
      
      r[:,i+1]=[y[0,i+1],y[1,i+1],y[2,i+1]]
      r_mag[i+1]=sqrt((r[0,i+1]**2)+(r[1,i+1]**2)+(r[2,i+1]**2))
      
      v[:,i+1]=[y[3,i+1],y[4,i+1],y[5,i+1]]
      
      acc[:,i+1]=-G*M*r[:,i+1]/(r_mag[i+1]**3)
      
      y_dot[:,i+1]=np.hstack([v[:,i+1],acc[:,i+1]])
      
      time[i+1]=time[i]+step
      
plt.plot(time,np.transpose(r))
      

