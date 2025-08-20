'''
The current version of the code has not been completed. The orbit is integrated by considering the J2 perturbations for
the calculation of acceleration. An objective function is created that calculates the 3D position difference for a slight 
change 'delta_y' to the initial state vector. Gradient descent method has been attempted to try and otpimize the orbit.
'''


from math import *
from sgp4.api import *
import numpy as np
import matplotlib.pyplot as plt
from sgp4.conveniences import sat_epoch_datetime
from datetime import datetime,timedelta
import time


mu=398600.4418E+9
R_e=6378.137E+3
G=6.667E-11
M=6E+24
J2=1.083E-3

def read_tle(filename):
#Function to read TLE file
    
    datContent = [i.strip() for i in open(filename).readlines()] 
    #Read TLE data
    
    s=datContent[1] 
    t=datContent[2]
    #Store TLE in separate variables to pass values to the relevant function
        
    satellite = Satrec.twoline2rv(s, t)
    #Convert TLE data to position and velocity
    
    start_epoch=satellite.jdsatepoch
    #Start epoch of the satellite is required when converting time_elapsed to Julian days
    
    return satellite, start_epoch


def propagate(filename, start, stop, step, output):
#Function to propagate orbit of the satellite within the given time_elapsed domain
    
    satellite, start_epoch = read_tle(filename)
    
    start_s=start_epoch*86400
    
    count=floor((stop-start)/step)
    #Number of iterations
    
    e=np.zeros(count)
    r=np.zeros([count,3])
    v=np.zeros([count,3])
    #Position and velocity domains
    
    time_elapsed=np.arange(start,stop,step)
    #time_elapsed domain 
    
    for i in range(count):
        
        jd=floor(start_epoch+time_elapsed[i]/86400)
        fr=(start_epoch+time_elapsed[i]/86400)%jd
        #Conversion of minutes to Julian days to pass values to the relevant sgp4 function
        
        e[i],a,b = satellite.sgp4(jd, fr)
        #Retrieving postion and velocity data for the corresponding Julian date in the time_elapsed domain
                    
        r[i,:]=np.asarray(a)*1000
        v[i,:]=np.asarray(b)*1000
        #Storing position and velocity in arrays to write in the output file
        
    out=np.column_stack((time_elapsed,r,v))    
    np.savetxt(output,out)
    #Writing data into the output file

def declare(count):
#Function to declare array variables
    
    time_elapsed=np.zeros(count)
    r=np.zeros([3,count])
    r_mag=np.zeros(count)
    v=np.zeros([3,count])
    acc=np.zeros([3,count])
    y=np.zeros([6,count])
    y_dot=np.zeros([6,count])
    diff=np.zeros([count,7])
    local_err=np.zeros(count)
    local_t_start=np.zeros(count)
    local_t_end=np.zeros(count)
    
    return time_elapsed, r, r_mag, v, acc, y, y_dot, diff, local_err, local_t_start, local_t_end

def acc_oblate(r):
#Function to calculate acceleration vector and magnitude with J2 perturbations
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    
    K1=(-3*mu*J2/2)*(R_e**2/r_mag**5)
    K2=r[2]**2/r_mag**2
    #Arbitrary constants to simplify final expression
    
    a_x=K1*r[0]*(1-(5*K2))
    a_y=K1*r[1]*(1-(5*K2))
    a_z=K1*r[2]*(3-(5*K2))
    
    acc=np.array([a_x,a_y,a_z])
    
    return acc

def acc_calc(r):
#Function to calculate acceleration vector and magnitude
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    acc=-G*M*r/(r_mag**3)
    
    return acc

def euler_int(filename, start, stop, step, delta_y, temp):
#Function to carry out integration using the Euler's method
    
    sgp4_data=np.loadtxt(filename,dtype=float)
    
    dim=np.shape(sgp4_data)
    count=dim[0]
    
    time_elapsed, r, r_mag, v, acc, y, y_dot, diff, local_err, local_t_start, local_t_end = declare(count)
    
    time_elapsed[0]=start

    x_init=sgp4_data[0,1]
    y_init=sgp4_data[0,2]
    z_init=sgp4_data[0,3]
    
    r[:,0]=[x_init,y_init,z_init]
    r_mag[0]=sqrt((x_init**2)+(y_init**2)+(z_init**2))
    #Initial position vector from sgp4 model
    
    x_dot_init=sgp4_data[0,4]
    y_dot_init=sgp4_data[0,5]
    z_dot_init=sgp4_data[0,6]
    v[:,0]=[x_dot_init,y_dot_init,z_dot_init]
    #Initial velocity vector from sgp4 model
    
    y[:,0]=np.hstack([r[:,0],v[:,0]])+delta_y
    
    acc[:,0]=-G*M*r[:,0]/(r_mag[0]**3)
    #Initial acceleration vector
    
    y_dot[:,0]=np.hstack([v[:,0],acc[:,0]])
    #Vectors used for iterations during integration
    
    begin=time.time()
    #Time value to be used for total time for integration
    
    diff_t=0
    #Variable to be used for storing average time per iteration
    
    for i in range(count-1):
        
        local_t_start[i]=time.time()
        #Variable to store start time of each iteration
        
        y[:,i+1]=y[:,i]+(step*y_dot[:,i])
      
        r[:,i+1]=[y[0,i+1],y[1,i+1],y[2,i+1]]
        r_mag[i+1]=sqrt((r[0,i+1]**2)+(r[1,i+1]**2)+(r[2,i+1]**2))
          
        v[:,i+1]=[y[3,i+1],y[4,i+1],y[5,i+1]]
          
        acc[:,i+1]=acc_calc(r[:,i+1])+acc_oblate(r[:,i+1])
          
        y_dot[:,i+1]=np.hstack([v[:,i+1],acc[:,i+1]])
          
        time_elapsed[i+1]=time_elapsed[i]+step
        
        local_t_end[i]=time.time()
        #Variable to store end time of each iteration
        
        diff_t=diff_t+local_t_end[i]-local_t_start[i]
        
    avg_t=diff_t/count
    end=time.time()
    
    output=np.column_stack((time_elapsed,np.transpose(r),np.transpose(v)))    
    
    print('\n time_elapsed per iteration: ',avg_t)
    print('\n Total time_elapsed: ',end-begin)
    
    return output, sgp4_data, count

def gradient_desc(filename, start, stop, step, delta_y):
#Function to optimize the orbit using gradient descent method

    no_it=100
    obj=np.zeros([no_it+1,6])
    temp=0
    
    for i in range(no_it):
        
        obj[i,:], sgp4_dat, integrator_dat = objective(filename, start, stop, step, delta_y, temp)
        
        grad=np.gradient(obj[i,:])
        delta_y=delta_y-(np.array([1e-9,1e-8,1e-8,1e-7,1e-6,1e-5])*step*grad) #cosntant decided based on trial and error
        temp=1
        
        if(abs(obj[i,0])<100): #Initial convergence criteria
            print('Converged!')
            break
        
        print('Objective: ',obj[i])
        
    plt.plot(sgp4_dat[:,0],sgp4_dat[:,2]-integrator_dat[:,2]) #Difference in Y position (Sample plot tested)
    #plt.plot(sgp4_dat[:,0],integrator_dat[:,2])
    
    return sgp4_dat, integrator_dat, delta_y, grad, obj
    
    
def objective(filename, start, stop, step, delta_y, temp):
#Function to calculate the objective function
    
    integrator_dat, sgp4_dat, count =euler_int(filename, start, stop, step, delta_y, temp)
    time_elapsed, r, r_mag, v, acc, y, y_dot, diff, local_err, local_t_start, local_t_end=declare(count)

    x_diff=sgp4_dat[count-1,1]-integrator_dat[count-1,1]
    y_diff=sgp4_dat[count-1,2]-integrator_dat[count-1,2]
    z_diff=sgp4_dat[count-1,3]-integrator_dat[count-1,3]
    x_dot_diff=sgp4_dat[count-1,4]-integrator_dat[count-1,4]
    y_dot_diff=sgp4_dat[count-1,5]-integrator_dat[count-1,5]
    z_dot_diff=sgp4_dat[count-1,6]-integrator_dat[count-1,6]
    
    obj=np.array([x_diff,y_diff,z_diff,x_dot_diff,y_dot_diff,z_dot_diff])
    #Objective function

    return obj, sgp4_dat, integrator_dat


def global_error_calc(filename, start, stop, step, y, temp):
#Function to calculate global error

    integrator_dat, sgp4_dat, count = euler_int(filename, start, stop, step, y, temp)
    
    diff=np.zeros([count,7])
    local_err=np.zeros(count)
    local_sum=0
    
    for i in range(count):
    
        diff[i,:]=sgp4_dat[i,:]-integrator_dat[i,:]
        local_err[i]=sqrt((diff[i,1]**2)+(diff[i,2]**2)+(diff[i,3]**2))
        local_sum=local_sum+(local_err[i]**2)
    
    global_err=sqrt((1/count)*local_sum)
    
    return global_err