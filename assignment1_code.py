'''
1. All variables written to an output file are in SI Units.
2. Variables ending with '_mag' imply magnitude of the variable.
3. Functions with more than one main plot need to be run separately to obtain each individual plot. While plotting, all but one plot command need to be commented.
'''


from math import *
from sgp4.api import *
import numpy as np
import matplotlib.pyplot as plt
from sgp4.conveniences import sat_epoch_datetime
from datetime import datetime,timedelta


mu=398600.4418E+9
R_e=6378.137E+3
omega_e=7292115e-11
f = 1/298.257223563

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
    #Start epoch of the satellite is required when converting time to Julian days
    
    return satellite, start_epoch


def time_conversion(filename, start):
#Function to convert time from Julian date to UTC, UT1 and GPS standard
    
    satellite, start_epoch = read_tle(filename)
    #Read TLE data by calling the appropriate function
    
    DUT_file = [i.strip() for i in open('finals2000A.daily').readlines()] 

    J2000=datetime(2000,1,1,12,0,0)
    #J2000 epoch 
    
    UTC = sat_epoch_datetime(satellite)
    UTC_time=UTC 
    #Variable that stores the datetime parameter
    
    UT1=UTC+timedelta(milliseconds=-19.6904)
    UT1_time=UT1
    
    GPS=UTC+timedelta(0,37)
    GPS_time=GPS
    
    UTC_J2000=(UTC.replace(tzinfo=None)-J2000).total_seconds()
    UTC_J2000_time=UTC_J2000+start
    #Start time of domain from the J2000 epoch in UTC
    
    UT1_J2000=(UT1.replace(tzinfo=None)-J2000).total_seconds()
    UT1_J2000_time=UT1_J2000+start
    #Start time of domain from the J2000 epoch in UT1
    
    GPS_J2000=(GPS.replace(tzinfo=None)-J2000).total_seconds()
    GPS_J2000_time=GPS_J2000+start
    #Start time of domain from the J2000 epoch in GPS
    
    return UTC_J2000


def propagate(filename, start, stop, step, output):
#Function to propagate orbit of the satellite within the given time domain
    
    satellite, start_epoch = read_tle(filename)
    UTC_J2000 = time_conversion(filename, start)
    #Function call to obtain satellite data, start epoch and time in terms of J2000
    
    count=floor((stop-start)/step)
    #Number of iterations
    
    e=np.zeros(count)
    r=np.zeros([count,3])
    v=np.zeros([count,3])
    #Position and velocity domains
    
    time=np.arange(start,stop,step)
    #Time domain 
    
    for i in range(count):
        
        jd=floor(start_epoch+time[i]/1440)
        fr=(start_epoch+time[i]/1440)%jd
        #Conversion of minutes to Julian days to pass values to the relevant sgp4 function
        
        e[i],a,b = satellite.sgp4(jd, fr)
        #Retrieving postion and velocity data for the corresponding Julian date in the time domain
                    
        r[i,:]=np.asarray(a)*1000
        v[i,:]=np.asarray(b)*1000
        #Storing position and velocity in arrays to write in the output file
        
        time[i]=time[i]+UTC_J2000
        #time in terms of J2000
        
    out=np.column_stack((time,r,v))    
    np.savetxt(output,out)
    #Writing data into the output file
    
    return count, output

    
def declare(count):
#Function to declare all relevant variables
    
    r_0=np.zeros(count)
    r_mag=np.zeros(count)
    v_mag=np.zeros(count)
    h=np.zeros([count,3])
    h_mag=np.zeros(count)
    e=np.zeros([count,3])
    e_mag=np.zeros(count)
    e_x=np.zeros([count,3])
    e_z=np.zeros([count,3])
    Phi=np.zeros(count)
    eta=np.zeros(count)
    a=np.zeros(count)
    inc=np.zeros(count)
    theta=np.zeros(count)
    r_dot_e=np.zeros(count)
    l1=np.zeros(count)
    l2=np.zeros(count)
    m1=np.zeros(count)
    m2=np.zeros(count)
    n1=np.zeros(count)
    n2=np.zeros(count)
    x=np.zeros(count)
    y=np.zeros(count)
    z=np.zeros(count)
    v_x=np.zeros(count)
    v_y=np.zeros(count)
    v_z=np.zeros(count)
    N_an=np.zeros([count,3])
    N_xy=np.zeros(count)
    N_an_mag=np.zeros(count)
    Omega=np.zeros(count)
    omega=np.zeros(count)
    r=np.zeros(count)
    r_xy=np.zeros(count)
    lambda_1=np.zeros(count)
    phi=np.zeros(count)
    lambda_0=np.zeros(count)
    
    return r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0


def cart_to_kep(filename, output):
#Function to convert cartesian to keplerian co-ordinates
    
    cart_data=np.loadtxt(filename, dtype=float)
    #Load cartesian co-ordinates

    dim=np.shape(cart_data)
    count=dim[0]
    
    r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0 = declare(count)

    time=cart_data[:,0]
    r=cart_data[:,[1,2,3]]
    v=cart_data[:,[4,5,6]]
    #Assigning cartesian co-ordinates to variables for the purpose of conversion

    for i in range(count):
     
        r_mag[i]=sqrt((r[i,0]**2)+(r[i,1]**2)+(r[i,2]**2))
        v_mag[i]=sqrt((v[i,0]**2)+(v[i,1]**2)+(v[i,2]**2))
        
        h[i,:]=np.cross(r[i,:],v[i,:])
        #Specific angular momentum vector
        
        h_mag[i]=sqrt((h[i,0]**2)+(h[i,1]**2)+(h[i,2]**2))
        
        a[i]=1/((2/r_mag[i])-((v_mag[i]**2)/mu))
        #Semi-major axis
        
        e[i,:]=(np.cross(v[i,:],h[i,:])/mu)-(r[i,:]/r_mag[i])
        #Eccentricity vector
        
        e_x[i,0]=e[i,0]
        e_z[i,2]=e[i,2]
        e_mag[i]=sqrt((e[i,0]**2)+(e[i,1]**2)+(e[i,2]**2))
        
        inc[i]=acos(h[i,2]/h_mag[i])*180/pi
        #Inclination angle
        
        r_dot_e[i]=np.dot((r[i,:]/r_mag[i]),e[i,:]/e_mag[i])
        #Dot product of position and eccentricty vector for calculating true anomaly
        
        if(np.dot((np.cross(e[i,:],r[i,:])),h[i,:])>0):
            sign_theta=1
        else:
            sign_theta=-1
        
        theta[i]=sign_theta*(acos(r_dot_e[i])*180/pi)
        #True anomaly with range of -180 degrees to +180 degrees
        
        N_an[i,:]=np.cross([0,0,1],h[i,:])
        #Ascending node vector
        
        N_xy[i]=sqrt((N_an[i,0]**2)+(N_an[i,1]))
        #Ascending node vector projection in the XY plane
        
        N_an_mag[i]=sqrt((N_an[i,0]**2)+(N_an[i,1]**2)+(N_an[i,2]**2))

        Omega[i]=atan2((N_an[i,1]/N_xy[i]),(N_an[i,0]/N_xy[i]))*180/pi
        #Right ascension
        
        omega[i]=acos(np.dot(N_an[i,:],e[i,:])/(N_an_mag[i]*e_mag[i]))*180/pi
        #Argument of perigee

    out=np.column_stack((time,a,e_mag,inc,Omega,omega,theta))    
    np.savetxt(output,out)
    
    #Start comment here
    plt.plot(time,theta)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('True anomaly (in degrees)', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,a)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('Semi-major axis (in metres)', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,inc)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('Inclination angle (in degrees)', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,e_mag)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('Eccentricity', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,omega)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('Argument of perigee (in degrees)', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,Omega)
    plt.xlabel('time (in minutes)', fontsize=16)
    plt.ylabel('Right ascension of ascending node (in degrees)', fontsize=16)
    plt.show()
    #End comment here

def kep_to_cart(filename, output, init_cart_file):
#Function to convert keplerian to cartesian co-ordinates

    mod_cart_file=output
    
    kep_data=np.loadtxt(filename, dtype=float)
    #Load keplerian co-ordinates
    
    dim=np.shape(kep_data)
    count=dim[0]
    
    r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0 = declare(count)
    
    time=kep_data[:,0]
    a=kep_data[:,1]
    e_mag=kep_data[:,2]
    inc=kep_data[:,3]*pi/180
    Omega=kep_data[:,4]*pi/180
    omega=kep_data[:,5]*pi/180
    theta=kep_data[:,6]*pi/180
    #Assigning keplerian co-ordinates to individual variables for conversion
    
    for i in range(count):
    
        r_mag[i]=a[i]*(1-e_mag[i]**2)/(1+(e_mag[i]*cos(theta[i])))
        
        phi[i]=r_mag[i]*cos(theta[i])
        #Cosine component of position vector
        
        eta[i]=r_mag[i]*sin(theta[i])
        #Sine component of position vector
    
        l1[i] = cos(Omega[i])*cos(omega[i])-sin(Omega[i])*sin(omega[i])*cos(inc[i])
        l2[i] = -cos(Omega[i])*sin(omega[i])-sin(Omega[i])*cos(omega[i])*cos(inc[i])
        m1[i] = sin(Omega[i])*cos(omega[i])+cos(Omega[i])*sin(omega[i])*cos(inc[i])
        m2[i] = -sin(Omega[i])*sin(omega[i])+cos(Omega[i])*cos(omega[i])*cos(inc[i])
        n1[i] = sin(omega[i])*sin(inc[i])
        n2[i] = cos(omega[i])*sin(inc[i])
        
        rot_mat1=[[l1[i],l2[i]],[m1[i],m2[i]],[n1[i],n2[i]]]
        rot_mat2=[[phi[i]],[eta[i]]]
        #Matrices used for conversion of co-ordinates
        
        [x[i],y[i],z[i]]=np.dot(rot_mat1,rot_mat2)
        #Cartesian position
        
        h_mag[i]=sqrt(mu*a[i]*(1-(e_mag[i]**2)))
        
        rot_mat3=[[-l1[i],l2[i]],[-m1[i],m2[i]],[-n1[i],n2[i]]]
        rot_mat4=[[sin(theta[i])],[e_mag[i]+cos(theta[i])]]
        
        [v_x[i],v_y[i],v_z[i]]=(mu/h_mag[i])*np.dot(rot_mat3,rot_mat4)
        #Cartesian Velocity
        
    out=np.column_stack((time,x,y,z,v_x,v_y,v_z))    
    np.savetxt(output,out)
       
    init_cart_dat1=np.loadtxt(init_cart_file)
    mod_cart_dat1=np.loadtxt(mod_cart_file)
    #Initial and modified cartesian data read for comparison
    
    diff_position=init_cart_dat1[:,[1,2,3]]-mod_cart_dat1[:,[1,2,3]]
    diff_velocity=init_cart_dat1[:,[4,5,6]]-mod_cart_dat1[:,[4,5,6]]
    #Residual error data
    
    fig, ax =plt.subplots(2)

    ax[0].plot(time,diff_position)
    ax[0].set_ylabel('Position (in metres)', fontsize=16)
    ax[0].legend(['X co-ordinate','Y co-ordinate','Z co-ordinate'], fontsize=16)
    
    ax[1].plot(time,diff_velocity)
    ax[1].set_xlabel('time (in minutes)', fontsize=16)
    ax[1].set_ylabel('Velocity (in m/s)', fontsize=16)
    ax[1].legend(['X velocity','Y velocity','Z velocity'], fontsize=16)
    plt.show()
    
def cart_to_spher(filename, output):
#Function to convert cartesian to spherical co-ordinates
    
    cart_data=np.loadtxt(filename, dtype=float)
    #Load cartesian co-ordinates
    
    dim=np.shape(cart_data)
    count=dim[0]
    
    r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0 = declare(count)

    time=cart_data[:,0]
    x=cart_data[:,1]
    y=cart_data[:,2]
    z=cart_data[:,3]
    #Assigning cartesian co-ordinates to individual variables for conversion
    
    R_p = (f-1)*R_e
    #Polar radius in terms of oblateness
    
    for i in range(count):
    
        r[i]=sqrt((x[i]**2)+(y[i]**2)+(z[i]**2))
        r_xy[i]=sqrt((x[i]**2)+(y[i]**2))
        
        lambda_1[i]=atan2((y[i]/r_xy[i]),(x[i]/r_xy[i]))*180/pi
        
        Phi[i]=asin(z[i]/r[i])*180/pi
        
        deno = (np.cos(lambda_1[i])**2)/R_e**2 + (np.sin(lambda_1[i])**2)/R_p**2
        r_0[i] = np.sqrt(1/deno)
        #Radius of Earth as a function of latitude
    
    #Start comment here
    plt.plot(time,Phi)
    plt.plot(time,lambda_1)
    plt.xlabel('Time (in minutes since J2000)', fontsize=16)
    plt.legend(['Latitude','Longitude'], fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(lambda_1,Phi)
    plt.ylabel('Latitude (in degrees)', fontsize=16)
    plt.xlabel('Longitude (in degrees)', fontsize=16)
    plt.show()
    #End comment here
    
    #Start comment here
    plt.plot(time,r-R_e)
    plt.plot(time,r-r_0)
    plt.legend(['For spherical Earth','For ellipsoid Earth'], fontsize=16)
    plt.xlabel('Time (in minutes since J2000)', fontsize=16)
    plt.ylabel('Altitude (in metres)', fontsize=16)
    plt.show()
    #End comment here
    
    out=np.column_stack((time,lambda_1,Phi,r))    
    np.savetxt(output,out)

def spher_to_cart(filename, output, init_cart_file):
#Function to convert spherical to cartesian co-ordinates 

    cart_data=np.loadtxt(filename, dtype=float)
    #Load spherical co-ordinates

    mod_cart_file=output
    
    dim=np.shape(cart_data)
    count=dim[0]
    
    r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0 = declare(count)
    
    time=cart_data[:,0]
    
    lambda_1=cart_data[:,1]*pi/180
    phi=cart_data[:,2]*pi/180
    r=cart_data[:,3]
    #Assigning spherical co-ordinates to individual variables for conversion
    
    for i in range(count):
        
        x[i]=r[i]*cos(phi[i])*cos(lambda_1[i])
        y[i]=r[i]*cos(phi[i])*sin(lambda_1[i])
        z[i]=r[i]*sin(phi[i])
            
    out=np.column_stack((time,x,y,z))    
    np.savetxt(output,out)
    
    init_cart_dat1=np.loadtxt(init_cart_file)
    mod_cart_dat2=np.loadtxt(mod_cart_file)
    #Storing initial and modified cartesian co-ordinates in separate variables for comparison
    
    diff_position=mod_cart_dat2[:,[1,2,3]]-init_cart_dat1[:,[1,2,3]]

    plt.plot(time,diff_position)
    plt.xlabel('Time (in minutes since J2000)', fontsize=16)
    plt.ylabel('Residual error in cartesian co-ordinates (in metres)', fontsize=16)
    plt.legend(['X co-ordinate','Y co-ordinate','Z co-ordinate'], fontsize=16)
    plt.show()
    
def co_rotating_spher(filename, output):
#Function to convert spherical to co-rotating spherical co-ordinates
    
    spher_data=np.loadtxt(filename, dtype=float)
    #Load inertial spherical co-ordinates
    
    dim=np.shape(spher_data)
    count=dim[0]
    
    r_0, r_mag, v_mag, h, h_mag, e, e_mag, e_x, e_z, Phi, eta, a, inc, theta, r_dot_e, m1, m2, l1, l2, n1, n2, x, y, z, v_x, v_y, v_z, N_an, N_xy, N_an_mag, Omega, omega, r, r_xy, lambda_1, phi, lambda_0 = declare(count)

    time=spher_data[:,0]
    step=time[1]-time[0]
    
    lambda_1=spher_data[:,1]*pi/180
    phi=spher_data[:,2]*pi/180
    r=spher_data[:,3]

    phi_AE=51.98988681616558
    lambda_AE=4.375752354503242
    #Location of Aerospace Faculty
     
    x_AE=R_e*cos(phi_AE)*cos(lambda_AE)
    y_AE=R_e*cos(phi_AE)*sin(lambda_AE)
    z_AE=R_e*sin(phi_AE)
    
    r_AE=sqrt((x_AE**2)+(y_AE**2)+(z_AE**2))
    #Cartesian co-ordinates of Aerospace Faculty to determine closest approach
    
    nearest_r=1e+10
    #Default large initial value 
    
    for i in range(count):
        
        lambda_0[i]=lambda_1[i]-(omega_e*step)
        
        temp=r[i]-r_AE
        #Temporary variable for distance between satellite and Faculty
        
        if(temp<nearest_r):

            nearest_r=temp
            near_time=time[i]
            near_time_MJD=((near_time-0.0196904)/864000)+2459907.5-2400000.5
        #To check if current distance between satellite and Faculty is the shortest
    
    out=np.column_stack((time,lambda_0,phi,r))    
    np.savetxt(output,out)

    plt.plot(lambda_0,phi)
    plt.plot(lambda_1,phi)
    plt.ylabel('Latitude (in degrees)', fontsize=16)
    plt.xlabel('Longitude (in degrees)', fontsize=16)
    plt.show()
    
    return near_time_MJD

