'''
~~~~~~~~~~~~~~~~LEGEND FOR choice~~~~~~~~~~~~~~~~~~~~~

1: No perturbations considered.
2. Perturbations due to J2 (oblateness) term considered.
3. Perturbations due to J2 and J2,2 terms considered.
4.Perutrbations due to atmospheric drag with a non co-rotating atmosphere considered.
5.Perutrbations due to atmospheric drag with a co-rotating atmosphere considered.
6. Perturbations due to solar radiation pressure considered with the region in the penumbra assumed to be a total eclipse.
7. Perturbations due to solar radiation pressure considered including the effect of solar radiation in the penumbra.
8. Perturbations due to the Sun and Moon are considered.


~~~~~~~~~~~~~~~~~LEGEND FOR plot_no~~~~~~~~~~~~~~~~~~~

1. Latitude v.s Longitude difference for perturbations due to J2
2. Acceleration due to J2 as a function of time
3. Latitude v.s Longitude difference for perturbations due to J2,2 compared to J2
4. Acceleration due to J2,2 as a function of time
5. Latitude v.s Longitude difference for perturbations due to atmospheric drag
6. Acceleration due to drag as a function of time
7. Latitude v.s Longitude difference for perturbations due to a co-rotating atmoshpere
8. Acceleration due to drag caused by a co-rotating atmosphere as a function of time
9. Latitude v.s Longitude difference for perturbations due to solar radiation pressure 
10. Acceleration due to solar radiation pressure as a function of time
11. Latitude v.s Longitude difference for perturbations due to solar radiation pressure (including effect of penumbra)
12. Acceleration due to solar radiation pressure (including penumbra) as a function of time
13. Latitude v.s Longitude difference for perturbations due to Sun and Moon
14. Acceleration due to Sun as a function of time
15. Acceleration due to Moon as a function of time

'''


from math import *
from sgp4.api import *
import numpy as np
import matplotlib.pyplot as plt

#List of constants used

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
    
    return count, output

def declare(count):
#Function to declare array variables
    
    time_elapsed=np.zeros(count)
    r=np.zeros([3,count])
    r_mag=np.zeros(count)
    v=np.zeros([3,count])
    acc=np.zeros([3,count])
    y=np.zeros([6,count])
    y_dot=np.zeros([6,count])
    diff_lambda_1=np.zeros([7,count])
    diff_Phi=np.zeros([7,count])
    diff_r=np.zeros([7,count])
    lambda_1=np.zeros(count)
    Phi=np.zeros(count)
    r_xy=np.zeros(count)
    
    return time_elapsed, r, r_mag, v, acc, y, y_dot, diff_lambda_1, diff_Phi, diff_r, lambda_1, Phi, r_xy
    
def acc_calc(r):
#Function to calculate acceleration vector and magnitude without perturbations
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    acc=-G*M*r/(r_mag**3)
    
    return acc

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

def acc_J22(r,lambda_1,Phi):
#Function to calculate acceleration vector and magnitude with J2,2 perturbations
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    K3=mu*J22*(R_e**2)/r_mag**4
    #Arbitrary constant
    
    a_r=-9*K3*(cos(Phi))**2*cos(2*(lambda_1-lambda_22))   
    a_phi=-3*K3*sin(2*Phi)*cos(2*(lambda_1-lambda_22))  
    a_lambda=-6*K3*cos(Phi)*sin(2*(lambda_1-lambda_22))
    
    acc_mag=sqrt((a_r**2)+(a_phi**2)+(a_lambda**2))
    a_x=acc_mag*cos(Phi)*cos(lambda_1)
    a_y=acc_mag*cos(Phi)*sin(lambda_1)
    a_z=acc_mag*sin(Phi)
    
    acc=np.array([a_x,a_y,a_z])
    
    return acc

def acc_drag(r,v):
#Function to calculate acceleration vector and magnitude with perturbations due to atmospheric drag without a co-rotating atmposphere
    
    h=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    v_mag=sqrt((v[0]**2)+(v[1]**2)+(v[2]**2))
    rho=rho_0*exp(-h/H)
    a_d=np.array(int(-0.5*rho*B*v_mag)*v)
    
    return a_d

def acc_drag_corotating(r,v0):
#Function to calculate acceleration vector and magnitude with perturbations due to atmospheric drag with a co-rotating atmposphere
    
    omega=np.array([0,0,omega_e])
    v=v0-np.cross(omega,r)
    
    h=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    v_mag=sqrt((v[0]**2)+(v[1]**2)+(v[2]**2))
    rho=rho_0*exp(-h/H)
    a_d=np.array(int(-0.5*rho*B*v_mag)*v)
    
    return a_d

def acc_SRP(r_es,choice):
#Function to calculate acceleration vector and magnitude with perturbations due to solar radiation pressure

    r_ss=r_es-r_sun
    #Subscript 'es' = Earth to Satellite, 'ss' = Sun to Satellite
    
    r_es_mag=sqrt((r_es[0]**2)+(r_es[1]**2)+(r_es[2]**2))
    r_ss_mag=sqrt((r_ss[0]**2)+(r_ss[1]**2)+(r_ss[2]**2))
    
    r_es_u=r_es/r_es_mag
    r_ss_u=r_ss/r_ss_mag
    
    r_ps=np.dot(r_ss_u,r_es)*r_ss_u
    r_ps_mag=sqrt((r_ps[0]**2)+(r_ps[1]**2)+(r_ps[2]**2))
    R_p=(r_ps_mag*R_sun)/r_ss_mag
    #Subscript 'ps' = Point P to Satellite
    
    r_ep=r_es-r_ps
    r_ep_mag=sqrt((r_ep[0]**2)+(r_ep[1]**2)+(r_ep[2]**2))
    h_g=r_ep_mag-R_e
    
    shadow_func=h_g/R_p
    
    if(choice==6):
        acc=umbra(shadow_func,r_ss_u)
    
    elif(choice==7):
        acc=penumbra(shadow_func,r_ss_u)
        
    return acc

def umbra(shadow_func,r_ss_u):
#Function to calculate acceleration vector and magnitude with perturbations due to solar radiation pressure without considering the effect of penumbra

    if(shadow_func<1):
        a_srp=0
    
    elif(shadow_func>=1):
        a_srp=-(B_r*W/c)*r_ss_u
        
    return a_srp

def penumbra(shadow_func,r_ss_u):
#Function to calculate acceleration vector and magnitude with perturbations due to solar radiation pressure considering the effect of penumbra

    if(shadow_func<-1):
        a_srp=0
    
    elif(shadow_func>=-1 and shadow_func<1):
        f_g=1-((1/pi)*acos(shadow_func))+((shadow_func/pi)*sqrt(1-shadow_func**2))
        f_a=1
        a_srp=-f_g*f_a*(B_r*W/c)*r_ss_u
        
    elif(shadow_func>=1):
        a_srp=-(B_r*W/c)*r_ss_u
    
    return a_srp
    

def acc_3BP(r,r_sun,r_moon):
#Function to calculate acceleration vector and magnitude with perturbations due to Sun and Moon
    
    r_mag=sqrt((r[0]**2)+(r[1]**2)+(r[2]**2))
    r_sun_mag=sqrt((r_sun[0]**2)+(r_sun[1]**2)+(r_sun[2]**2))
    r_moon_mag=sqrt((r_moon[0]**2)+(r_moon[1]**2)+(r_moon[2]**2))
    
    a_x_sun=mu_sun*(((r_sun[0]-r[0])/(r_sun_mag-r_mag)**3)-(r_sun[0]/r_sun_mag**3)) 
    a_y_sun=mu_sun*(((r_sun[1]-r[1])/(r_sun_mag-r_mag)**3)-(r_sun[1]/r_sun_mag**3))
    a_z_sun=mu_sun*(((r_sun[2]-r[2])/(r_sun_mag-r_mag)**3)-(r_sun[2]/r_sun_mag**3)) 
    
    a_x_moon=mu_moon*(((r_moon[0]-r[0])/(r_moon_mag-r_mag)**3)-(r_moon[0]/r_moon_mag**3))
    a_y_moon=mu_moon*(((r_moon[1]-r[1])/(r_moon_mag-r_mag)**3)-(r_moon[1]/r_moon_mag**3))
    a_z_moon=mu_moon*(((r_moon[2]-r[2])/(r_moon_mag-r_mag)**3)-(r_moon[2]/r_moon_mag**3))
    
    acc_sun=np.array([a_x_sun,a_y_sun,a_z_sun])
    acc_moon=np.array([a_x_moon,a_y_moon,a_z_moon])
    
    return acc_sun, acc_moon
    #Returning acceleration of Sun and Moon for separate plots

def euler_int(filename, start, stop, step, choice):
#Function to carry out integration using the Euler's method based input variable 'choice' to account for different perturbations
    
    sgp4_data=np.loadtxt(filename, dtype=float)
    
    dim=np.shape(sgp4_data)
    count=dim[0]
    
    time_elapsed, r, r_mag, v, acc, y, y_dot, diff_lambda_1, diff_Phi, diff_r, lambda_1, Phi, r_xy = declare(count)

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
    
    acc[:,0]=-G*M*r[:,0]/(r_mag[0]**3)
    #Initial acceleration vector
    
    y[:,0]=np.hstack([r[:,0],v[:,0]])
    y_dot[:,0]=np.hstack([v[:,0],acc[:,0]])
    #Vectors used for iterations during integration
    
    r_xy[0]=sqrt((x_init**2)+(y_init**2))
    lambda_1[0]=atan2((y_init/r_xy[0]),(x_init/r_xy[0]))*180/pi        
    Phi[0]=asin(z_init/r_mag[0])*180/pi
    
    acc_s=np.zeros(count)
    acc_m=np.zeros(count)
        
    for i in range(count-1):
        
        y[:,i+1]=y[:,i]+(step*y_dot[:,i])
      
        r[:,i+1]=[y[0,i+1],y[1,i+1],y[2,i+1]]
        r_mag[i+1]=sqrt((r[0,i+1]**2)+(r[1,i+1]**2)+(r[2,i+1]**2))
          
        v[:,i+1]=[y[3,i+1],y[4,i+1],y[5,i+1]]
        
        r_xy[i+1]=sqrt((r[0,i+1]**2)+(y[1,i+1]**2))
        
        Atan=atan2((r[1,i+1]/r_xy[i+1]),(r[0,i+1]/r_xy[i+1]))
        
        lambda_1[i+1]=Atan*180/pi
        
        Phi[i+1]=asin(r[2,i+1]/r_mag[i+1])*180/pi
          
        if(choice==1):
            acc[:,i+1]=acc_calc(r[:,i+1])
            
        elif(choice==2):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_oblate(r[:,i+1])
          
        elif(choice==3):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_oblate(r[:,i+1])+acc_J22(r[:,i+1],lambda_1[i+1],Phi[i+1])
            
        elif(choice==4):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_drag(r[:,i+1],v[:,i+1])
            
        elif(choice==5):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_drag_corotating(r[:,i+1],v[:,i+1])
            
        elif(choice==6):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_SRP(r[:,i+1],choice)
            
        elif(choice==7):
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_SRP(r[:,i+1],choice)
             
        elif(choice==8):            
            acc_sun, acc_moon = acc_3BP(r[:,i+1], r_sun, r_moon)
            acc[:,i+1]=acc_calc(r[:,i+1])+acc_sun+acc_moon
            
            acc_sun_mag=sqrt((acc_sun[0]**2)+(acc_sun[1]**2)+(acc_sun[2]**2))
            acc_s[i+1]=acc_sun_mag
            
            acc_moon_mag=sqrt((acc_moon[0]**2)+(acc_moon[1]**2)+(acc_moon[2]**2))
            acc_m[i+1]=acc_moon_mag
            
            if(i==0):
                acc_s[i]=acc_s[i+1]
                acc_m[i]=acc_m[i+1]
            
        y_dot[:,i+1]=np.hstack([v[:,i+1],acc[:,i+1]])
        
        time_elapsed[i+1]=time_elapsed[i]+step
         
    output=np.column_stack((time_elapsed,lambda_1,Phi,np.transpose(r)))  
    to_write=np.column_stack((time_elapsed,np.transpose(r),np.transpose(v)))
    
    if(choice==8):
        return output,to_write,acc,acc_s,acc_m
        #Passing acceleration of Sun and Moon for plotting
    
    else:
        return output,to_write,acc


def plot(filename, start, stop, step, plot_no):
#Function to plot variousparameters based on the input variable 'plot_no'
    
    output_1, file_1, acc_1 = euler_int(filename, start, stop, step, 1)
    output_2, file_2, acc_2 = euler_int(filename, start, stop, step, 2)
    output_3, file_3, acc_3 = euler_int(filename, start, stop, step, 3)
    output_4, file_4, acc_4 = euler_int(filename, start, stop, step, 4)
    output_5, file_5, acc_5 = euler_int(filename, start, stop, step, 5)
    output_6, file_6, acc_6 = euler_int(filename, start, stop, step, 6)
    output_7, file_7, acc_7 = euler_int(filename, start, stop, step, 7)
    output_8, file_8, acc_8, acc_s, acc_m = euler_int(filename, start, stop, step, 8)
    
    if(plot_no==0):
        write_dat(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8)
        #To write output files without plotting
    
    dim=np.shape(output_1)
    count=dim[0]
    
    acc_1_mag=np.zeros(count)
    acc_2_mag=np.zeros(count)
    acc_3_mag=np.zeros(count)
    acc_4_mag=np.zeros(count)
    acc_5_mag=np.zeros(count)
    acc_6_mag=np.zeros(count)
    acc_7_mag=np.zeros(count)
    acc_8_mag=np.zeros(count)
    
    time_elapsed, r, r_mag, v, acc, y, y_dot, diff_lambda_1, diff_Phi, diff_r, lambda_1, Phi, r_xy = declare(count)
    
    diff_r[0,:]=output_1[:,3]-output_2[:,3]
    diff_lambda_1[0,:]=output_1[:,1]-output_2[:,1]
    diff_Phi[0,:]=output_1[:,2]-output_2[:,2]
    time_elapsed=output_1[:,0]
    
    diff_r[1,:]=output_3[:,3]-output_2[:,3]
    diff_lambda_1[1,:]=output_3[:,1]-output_2[:,1]
    diff_Phi[1,:]=output_3[:,2]-output_2[:,2]
    
    diff_r[2,:]=output_1[:,3]-output_4[:,3]
    diff_lambda_1[2,:]=output_1[:,1]-output_4[:,1]
    diff_Phi[2,:]=output_1[:,2]-output_4[:,2]
    
    diff_r[3,:]=output_5[:,3]-output_4[:,3]
    diff_lambda_1[3,:]=output_5[:,1]-output_4[:,1]
    diff_Phi[3,:]=output_5[:,2]-output_4[:,2]
    
    diff_r[4,:]=output_1[:,3]-output_6[:,3]
    diff_lambda_1[4,:]=output_1[:,1]-output_6[:,1]
    diff_Phi[4,:]=output_1[:,2]-output_6[:,2]
    
    diff_r[5,:]=output_7[:,3]-output_6[:,3]
    diff_lambda_1[5,:]=output_7[:,1]-output_6[:,1]
    diff_Phi[5,:]=output_7[:,2]-output_6[:,2]
        
    diff_r[6,:]=output_1[:,3]-output_8[:,3]
    diff_lambda_1[6,:]=output_1[:,1]-output_8[:,1]
    diff_Phi[6,:]=output_1[:,2]-output_8[:,2]
    
    for i in range(count):
        
        acc_1_mag[i]=sqrt((acc_1[0,i]**2)+(acc_1[1,i]**2)+(acc_1[2,i]**2))
        acc_2_mag[i]=sqrt((acc_2[0,i]**2)+(acc_2[1,i]**2)+(acc_2[2,i]**2))
        acc_3_mag[i]=sqrt((acc_3[0,i]**2)+(acc_3[1,i]**2)+(acc_3[2,i]**2))
        acc_4_mag[i]=sqrt((acc_4[0,i]**2)+(acc_4[1,i]**2)+(acc_4[2,i]**2))
        acc_5_mag[i]=sqrt((acc_5[0,i]**2)+(acc_5[1,i]**2)+(acc_5[2,i]**2))
        acc_6_mag[i]=sqrt((acc_6[0,i]**2)+(acc_6[1,i]**2)+(acc_6[2,i]**2))
        acc_7_mag[i]=sqrt((acc_7[0,i]**2)+(acc_7[1,i]**2)+(acc_7[2,i]**2))
        acc_8_mag[i]=sqrt((acc_8[0,i]**2)+(acc_8[1,i]**2)+(acc_8[2,i]**2))
        
    if(plot_no==1):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to J2', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[0,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.set_ylim([-10,10])
        ax.plot(time_elapsed,diff_Phi[0,:],color='red')
        ax2.plot(time_elapsed,diff_r[0,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==2):
        plt.plot(time_elapsed,acc_2_mag-acc_1_mag)
        plt.title('Acceleration due to J2 as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to J2 perturbations (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==3):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to J2,2 compared to J2', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[1,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.set_ylim([-10,10])
        ax.plot(time_elapsed,diff_Phi[1,:],color='red')
        ax2.plot(time_elapsed,diff_r[1,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==4):
        plt.plot(time_elapsed,acc_3_mag-acc_2_mag)
        plt.title('Acceleration due to J2,2 as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to J2,2 perturbations (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==5):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to atmospheric drag', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[2,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.plot(time_elapsed,diff_Phi[2,:],color='red')
        ax2.plot(time_elapsed,diff_r[2,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==6):
        plt.plot(time_elapsed,acc_4_mag-acc_1_mag)
        plt.title('Acceleration due to drag as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to drag (m/s^2)', fontsize=16)
        plt.show()
    
    elif(plot_no==7):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to a co-rotating atmoshpere', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[3,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.plot(time_elapsed,diff_Phi[3,:],color='red')
        ax2.plot(time_elapsed,diff_r[3,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==8):
        plt.plot(time_elapsed,acc_5_mag-acc_4_mag)
        plt.title('Acceleration due to a co-rotating atmoshpere as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to co-rotating atmospheric drag (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==9):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to solar radiation pressure', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[4,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.plot(time_elapsed,diff_Phi[4,:],color='red')
        ax2.plot(time_elapsed,diff_r[4,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==10):
        plt.plot(time_elapsed,acc_6_mag-acc_1_mag)
        plt.title('Acceleration due to solar radiation pressure as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to Solar Radiation Pressure (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==11):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to solar radiation pressure (including effect of penumbra)', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[5,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.plot(time_elapsed,diff_Phi[5,:],color='red')
        ax2.plot(time_elapsed,diff_r[5,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==12):
        plt.plot(time_elapsed,acc_7_mag-acc_6_mag)
        plt.title('Acceleration due to solar radiation pressure (including effect of penumbra) as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to Solar Radiation Pressure including penumbra (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==13):
        fig, ax=plt.subplots()
        ax.set_title('Latitude v.s Longitude difference for perturbations due to Sun and Moon', fontsize=16)
        ax2=ax.twinx()
        ax.plot(time_elapsed,diff_lambda_1[6,:])
        ax.set_xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        ax.set_ylabel('Latitude/Longitude difference (in degrees)', fontsize=16)
        ax.plot(time_elapsed,diff_Phi[6,:],color='red')
        ax2.plot(time_elapsed,diff_r[6,:],color='black')
        ax2.set_ylabel('Radial difference (in metres)', fontsize=16)
        fig.legend(['latitude','longitude','radial'], fontsize=16)
        plt.show()
        
    elif(plot_no==14):
        plt.plot(time_elapsed,acc_s)
        plt.title('Acceleration due to the Sun as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to perturbations caused by Sun (m/s^2)', fontsize=16)
        plt.show()
        
    elif(plot_no==15):
        plt.plot(time_elapsed,acc_m)
        plt.title('Acceleration due to the moon as a function of time', fontsize=16)
        plt.xlabel('Time elapsed (in seconds since TLE epoch)', fontsize=16)
        plt.ylabel('Acceleration due to perturbations caused by Moon (m/s^2)', fontsize=16)
        plt.show()
        

def write_dat(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8):
#Function to write the output into .dat files
    
    np.savetxt('output_1.dat',file_1)
    np.savetxt('output_2.dat',file_2)
    np.savetxt('output_3.dat',file_3)
    np.savetxt('output_4.dat',file_4)
    np.savetxt('output_5.dat',file_5)
    np.savetxt('output_6.dat',file_6)
    np.savetxt('output_7.dat',file_7)
    np.savetxt('output_8.dat',file_8)

