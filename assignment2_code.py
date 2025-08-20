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
    #Start epoch of the satellite is required when converting time_elpased to Julian days
    
    return satellite, start_epoch


def propagate(filename, start, stop, step, output):
#Function to propagate orbit of the satellite within the given time_elpased domain
    
    satellite, start_epoch = read_tle(filename)
    
    start_s=start_epoch*86400
    
    count=floor((stop-start)/step)
    #Number of iterations
    
    e=np.zeros(count)
    r=np.zeros([count,3])
    v=np.zeros([count,3])
    #Position and velocity domains
    
    time_elpased=np.arange(start,stop,step)
    #time_elpased domain 
    
    for i in range(count):
        
        jd=floor(start_epoch+time_elpased[i]/86400)
        fr=(start_epoch+time_elpased[i]/86400)%jd
        #Conversion of minutes to Julian days to pass values to the relevant sgp4 function
        
        e[i],a,b = satellite.sgp4(jd, fr)
        #Retrieving postion and velocity data for the corresponding Julian date in the time_elpased domain
                    
        r[i,:]=np.asarray(a)*1000
        v[i,:]=np.asarray(b)*1000
        #Storing position and velocity in arrays to write in the output file
        
    out=np.column_stack((time_elpased,r,v))    
    np.savetxt(output,out)
    #Writing data into the output file
    
    return count, output

def declare(count):
#Function to declare array variables
    
    time_elpased=np.zeros(count)
    r=np.zeros([3,count])
    r_mag=np.zeros(count)
    v=np.zeros([3,count])
    acc=np.zeros([3,count])
    y=np.zeros([6,count])
    y_dot=np.zeros([6,count])
    diff_1=np.zeros([count,7])
    diff_2=np.zeros(7)
    local_err_1=np.zeros(count)
    local_err_2=np.zeros(count)
    local_t_start=np.zeros(count)
    local_t_end=np.zeros(count)
    
    return time_elpased, r, r_mag, v, acc, y, y_dot, diff_1, diff_2, local_err_1, local_err_2, local_t_start, local_t_end
    
def acc_calc(r):
#Function to calculate acceleration vector and magnitude
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    acc=-G*M*r/(r_mag**3)
    
    return acc

def euler_int(filename, start, stop, step):
#Function to carry out integration using the Euler's method
    
    sgp4_data=np.loadtxt(filename, dtype=float)
    
    dim=np.shape(sgp4_data)
    count=dim[0]
    
    time_elpased, r, r_mag, v, acc, y, y_dot, diff_1, diff_2, local_err_1, local_err_2, local_t_start, local_t_end = declare(count)
    
    time_elpased[0]=start

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
    
    acc[:,0]=-G*M*r[:,0]/(r_mag[0]**3)
    #Initial acceleration vector
    
    y[:,0]=np.hstack([r[:,0],v[:,0]])
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
          
        acc[:,i+1]=acc_calc(r[:,i+1])
          
        y_dot[:,i+1]=np.hstack([v[:,i+1],acc[:,i+1]])
          
        time_elpased[i+1]=time_elpased[i]+step
        
        local_t_end[i]=time.time()
        #Variable to store end time of each iteration
        
        diff_t=diff_t+local_t_end[i]-local_t_start[i]
        
    avg_t=diff_t/count
    end=time.time()
    
    output=np.column_stack((time_elpased,np.transpose(r),np.transpose(v)))    
    
    print('\n time_elpased per iteration: ',avg_t)
    print('\n Total time_elpased: ',end-begin)
    
    return output


def adams_bashforth(filename, start, stop, step):
#Function used to carry out integration using Adams Bashforth method
    
    sgp4_data=np.loadtxt(filename, dtype=float)
    euler_dat=euler_int(filename, start, stop, step)
    #Euler's method used for initial 'm' values
    
    dim=np.shape(sgp4_data)
    count=dim[0]
    
    time_elpased, r, r_mag, v, acc, y, y_dot, diff_1, diff_2, local_err_1, local_err_2, local_t_start, local_t_end = declare(count)
    
    time_elpased[0]=start

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
    
    acc[:,0]=-G*M*r[:,0]/(r_mag[0]**3)
    #Initial acceleration vector
    
    for k in range(5):
        
        acc_data=acc_calc(np.array([euler_dat[k,1],euler_dat[k,2],euler_dat[k,3]]))
        y[:,k]=[euler_dat[k,1],euler_dat[k,2],euler_dat[k,3],euler_dat[k,4],euler_dat[k,5],euler_dat[k,6]]
        y_dot[:,k]=[euler_dat[k,4],euler_dat[k,5],euler_dat[k,6], acc_data[0], acc_data[1], acc_data[2]]
        #Initial values obtained from Euler's method
        
    for i in range(count):
        
        if(i<=4):
            continue
        
        phi_AB2=0.5*((3*y_dot[:,i])-y_dot[:,i-1])
        phi_AB3=(1/12)*((23*y_dot[:,i])-(16*y_dot[:,i-1])+(5*y_dot[:,i-2]))
        phi_AB4=(1/24)*((55*y_dot[:,i])-(59*y_dot[:,i-1])+(37*y_dot[:,i-2])-(9*y_dot[:,i-3]))
        phi_AB5=(1/710)*((1901*y_dot[:,i])-(2774*y_dot[:,i-1])+(2616*y_dot[:,i-2])-(1274*y_dot[:,i-3])+(251*y_dot[:,i-4]))
        
        if(i<count-1):
            
            y[:,i+1]=y[:,i]+(step*phi_AB5)
            acc_data=acc_calc(np.array([y[0,i+1],y[1,i+1],y[2,i+1]]))
            y_dot[:,i+1]=[y[3,i+1],y[4,i+1],y[5,i+1],acc_data[0],acc_data[1],acc_data[2]]

        
def local_error_calc(filename, start, stop, step_1, step_2, sgp4_1, sgp4_2):
#Function to calculate the local error and compare for two step sizes
    
    count_1, sgp4_1 = propagate(filename, start, stop ,step_1, sgp4_1)
    count_2, sgp4_2 = propagate(filename, start, stop ,step_2, sgp4_2)
    #Sgp4 model needs to be run with the same step size beforehand
    
    integrator_1_dat = euler_int(sgp4_1, start, stop, step_1)
    integrator_2_dat = euler_int(sgp4_2, start, stop, step_2)
    
    sgp4_1_dat=np.loadtxt(sgp4_1, dtype=float)
    sgp4_2_dat=np.loadtxt(sgp4_2, dtype=float)
    
    time_elpased_1, r, r_mag, v, acc, y, y_dot, diff_1, diff_2, local_err_1, local_err_2, local_t_start, local_t_end = declare(count_1)
    
    time_elpased_1=sgp4_1_dat[:,0]
    
    for i in range(count_1):
        
        diff_1[i,:]=sgp4_1_dat[i,:]-integrator_1_dat[i,:]
        local_err_1[i]=sqrt((diff_1[i,1]**2)+(diff_1[i,2]**2)+(diff_1[i,3]**2)+(diff_1[i,4]**2)+(diff_1[i,5]**2)+(diff_1[i,6]**2))
        #local error expressed as rms value of position and velocity difference
        
    plt.plot(time_elpased_1,local_err_1)
    
    time_elpased_2, r, r_mag, v, acc, y, y_dot, diff_1, diff_2, local_err_1, local_err_2, local_t_start, local_t_end = declare(count_2)
    
    time_elpased_2=sgp4_2_dat[:,0]
    
    for i in range(count_2):
        
        diff_2=sgp4_2_dat[i,:]-integrator_2_dat[i,:]
        local_err_2[i]=sqrt((diff_2[1]**2)+(diff_2[2]**2)+(diff_2[3]**2)+(diff_2[4]**2)+(diff_2[5]**2)+(diff_2[6]**2))
    
    plt.plot(time_elpased_2,local_err_2)
    plt.legend(['Step size 1','Step size 2'])
    plt.xlabel('Time elapsed in seconds since the TLE epoch', fontsize=16)
    plt.ylabel('Local error', fontsize=16)
    plt.show()
    
def global_error_calc(filename, start, stop, step, sgp4_out, integrator_out):
#Function to calculate global error and plot the effect of step size on it

    count=np.zeros(len(step))
    global_err=np.zeros(len(step))
    
    for i in range(len(step)):
        
        count[i], sgp4 = propagate(filename, start, stop ,step[i], sgp4_out[i])
        integrator_dat = euler_int(sgp4_out[i], start, stop, step[i])
        
        sgp4_dat=np.loadtxt(sgp4, dtype=float)
        
        diff=np.zeros([int(count[i]),7])
        local_err=np.zeros(int(count[i]))
        local_sum=0
        
        for j in range(int(count[i])):
        
            diff[j,:]=sgp4_dat[j,:]-integrator_dat[j,:]
            local_err[j]=sqrt((diff[j,1]**2)+(diff[j,2]**2)+(diff[j,3]**2)+(diff[j,4]**2)+(diff[j,5]**2)+(diff[j,6]**2))
            local_sum=local_sum+(local_err[j]**2)
        
        global_err[i]=sqrt((1/count[i])*local_sum)
    
    plt.plot(step,global_err)
    plt.xlabel('Step size (in seconds)', fontsize=16)
    plt.ylabel('Global RMS Error', fontsize=16)
    plt.show()
    
def local_stability(filename, tweak):
#Function to calculate and plot 3D local position difference
    
    sgp4_data=np.loadtxt(filename, dtype=float)
    
    dim=np.shape(sgp4_data)
    count=dim[0]
    
    diff=np.zeros([count,7])
    init=np.zeros(7)
    local_err=np.zeros([count,6])
    
    time_elapsed=sgp4_data[:,0]
    start=sgp4_data[0,0]
    stop=sgp4_data[count-1,0]
    step=sgp4_data[1,0]-sgp4_data[0,0]
    
    integrator_init_dat=euler_int(filename, start, stop, step)

    for i in range(6):
        
        init=sgp4_data[0,:]
        init[i+1]=init[i+1]*(1+tweak[i]/100)
        #Tweaked orbit value calculated using % value as input for 'tweak' variable
        
        sgp4_data[0,:]=init
        
        np.savetxt('temp_out.dat',sgp4_data)
        
        integrator_dat=euler_int('temp_out.dat', start, stop, step)
                
        for j in range(count):
                
                diff[j,:]=integrator_dat[j,:]-integrator_init_dat[j,:]
                local_err[j,i]=sqrt((diff[j,1]**2)+(diff[j,2]**2)+(diff[j,3]**2))
        
        #start comment here while running global_stability function    
        plt.plot(time_elapsed,local_err[:,i])        
        plt.xlabel('Time in seconds since TLE epoch', fontsize=16)
        plt.ylabel('Local 3D Position Difference (in metres)', fontsize=16)
        plt.legend(['x-Position','y-Position','z-Position','x-Velocity','y-Velocity','z-Velocity'])
 
    plt.show()
    #end comment here
             
    return local_err
              
     
def global_stability(filename, no_tweaks):
#Function to calculate and plot 3D global position difference
    
        local_sum_x=np.zeros(no_tweaks); local_sum_y=np.zeros(no_tweaks); local_sum_z=np.zeros(no_tweaks)
        local_sum_vx=np.zeros(no_tweaks); local_sum_vy=np.zeros(no_tweaks); local_sum_vz=np.zeros(no_tweaks)
        gloabl_rms_x=np.zeros(no_tweaks); gloabl_rms_y=np.zeros(no_tweaks); gloabl_rms_z=np.zeros(no_tweaks)
        gloabl_rms_vx=np.zeros(no_tweaks); gloabl_rms_vy=np.zeros(no_tweaks); gloabl_rms_vz=np.zeros(no_tweaks);
        
        tweak_count=np.zeros(no_tweaks)
        
        tweak=np.zeros([no_tweaks,6])
        #Arbitrary number of tweaks provided (no_tweaks) because of no specific requirement
        
        for i in range(no_tweaks):
            
            tweak_count[i]=i+1
            tweak[i,:]=[i,i,i,i,i,i]
            #Tweaks to vary from 1% to 10% of the original value(s)
            
            local_err=local_stability(filename, tweak[i,:])
            
            dim=np.shape(local_err)
            count=dim[0]
    
            local_sum_x[i]=sum(local_err[:,0]**2)
            local_sum_y[i]=sum(local_err[:,1]**2)
            local_sum_z[i]=sum(local_err[:,2]**2)
            local_sum_vx[i]=sum(local_err[:,3]**2)
            local_sum_vy[i]=sum(local_err[:,4]**2)
            local_sum_vz[i]=sum(local_err[:,5]**2)
            
            gloabl_rms_x[i]=(1/count)*sqrt(local_sum_x[i])
            gloabl_rms_y[i]=(1/count)*sqrt(local_sum_y[i])
            gloabl_rms_z[i]=(1/count)*sqrt(local_sum_z[i])
            gloabl_rms_vx[i]=(1/count)*sqrt(local_sum_vx[i])
            gloabl_rms_vy[i]=(1/count)*sqrt(local_sum_vy[i])
            gloabl_rms_vz[i]=(1/count)*sqrt(local_sum_vz[i])
            
        plt.plot(tweak_count, gloabl_rms_x)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()
        
        plt.plot(tweak_count, gloabl_rms_y)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()
        
        plt.plot(tweak_count, gloabl_rms_z)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()
        
        plt.plot(tweak_count, gloabl_rms_vx)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()
        
        plt.plot(tweak_count, gloabl_rms_vy)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()
        
        plt.plot(tweak_count, gloabl_rms_vz)
        plt.xlabel('% Tweak in orbit', fontsize=16)
        plt.ylabel('Global 3D position difference (in metres)', fontsize=16)
        plt.show()