# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:22:25 2022

@author: jain_
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
In this program, we are going to look at how to calculate the 
relationship between "stuck time" and the "potential depth" in a 
non-equilibrium situation. Here, a 2 micrometer sized particle is stuck in a 
circular trap of radius ~8um. On this trap, there is a tangential drive force 
and a cosine potential of depth of the order kbT created by an optical tweezer 
and inside water medium. The drive force is of the order of ~1 pN. 
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# plt.style.use('ggplot')
# import sympy as sp
import math
import time
#%%
F0 = 400e-15# Drive force
T=300         #T #Kelvin is 25degreeCelsius
eta=8.9e-4    #eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
pi=np.pi;    
a=2e-6        # Micrometers Size of particle
r=10e-6        #Radius of the outer circle
periodR = 2*pi*r;
k_b=1.38e-23
zeta=3*6*pi*eta*a  #A factor 3 is needed to account for the plate interaction
m=4.0*pi/3*a**3*1100 #density of polystyrene is 1100kg/m^3
kBT = k_b*T;
#%
#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
V_p=20-6; #micrometers per second
# Initial Conditions: Potential well

rMat= np.linspace(0.0,periodR,1000);
thetaMat = np.linspace(0.0,2*pi,1000)
#Potential Function
Ufunc_ = 1000*k_b*T*np.cos(thetaMat) 
Ufunc_ = Ufunc_- max(Ufunc_)
Uforce_ = -np.diff(Ufunc_)/np.diff(rMat)


fig, axs = plt.subplots(3)
fig.suptitle('Potential, force, and tilted potential')
axs[0].plot(rMat,Ufunc_)

rMat_ = np.linspace(0.0,periodR,999)
axs[1].plot(rMat_,Uforce_)

from scipy.interpolate import interp1d

Ufunc = interp1d(rMat,Ufunc_)
Uforce = interp1d(rMat_,Uforce_)

axs[2].plot(rMat/periodR*2*pi,Ufunc_-F0*rMat)
# We calculate deltaU and the minima position
#%% Calculation
nums=int(1e5) # Simulation steps
dt =0.001     # Delta t of the simulation (Could be variable)
nump=4     # number of particles

dim  = 1   # system dimension (x,y,z)
std  =  np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
# np.random.seed(7) # initialize random number generator with a seed=0
R = np.zeros([nump,dim]) # array for starting & current positions  
  
for ii in range(nump):
    R[ii][0]=periodR*ii/2/nump;

#for i in range(nump):
#    R[i,0] = 3.72e-5; #i/nump*periodR

W = np.zeros([nump,dim]) # array to store current random forcces
F = np.zeros([nump,dim]) # array to store external force
Rs = np.zeros([nums,nump,dim]) # array to store positions at all steps
#Ws = np.zeros([nums,nump,dim]) # array to store random forces at all steps
Fs = np.zeros([nums,nump,dim]) # array to store external forces at all steps

timeMat  = np.zeros([nums]) # an array to store time at all steps
#%%
for i in range(nums): # repeat the following operations from i=0 to nums-1
    if (R.any()>periodR or R.any()<0):
        pR = periodR*1e20;
        R__=((R*1e20)%pR)
        R_=R__/1e20;
        F = Uforce(R_) +F0#*np.ones([1,nump])
    else:
        F = Uforce(R) +F0#*np.ones([1,nump])
        
        
    
    W = std*np.random.randn(nump,dim) # generate an array of random forces

    F+=W/dt; # Added random noise force to trap force


    
    # Now we let the force transmit like in Newton's cradle
    F_new=F;
    for kk in range(nump):  # We have to let it transmit through the chain of nump particles
        for ii in range(nump): #If it can't jump, it simply transmits force
    
            if(abs(R[ii][0]%periodR-R[ii-1][0]%periodR)<1.05*a):

                if(R[ii][0]%periodR>R[ii-1][0]%periodR):
                    F[ii][0]=+2*F0
                    F[ii-1][0]=-2*F0
                else:
                    F[ii-1][0]=-2*F0
                    F[ii][0]=+2*F0
                       
    
            if(R[ii][0]%periodR-R[ii-1][0]%periodR>0 and abs(R[ii][0]%periodR-R[ii-1][0]%periodR-2*a)<(F[ii-1][0]-F[ii][0])*dt/zeta):
                    F_new[ii-1][0]=0
                    F_new[ii][0] += F[ii-1][0]    
                    
            if(R[ii][0]%periodR-R[ii-1][0]%periodR>0 and F[ii][0]<0):
                if F[ii-1][0]<0 and abs(F[ii][0]-F[ii-1][0])*dt/zeta>2*a :
                    F_new[ii][0]=0
                    F_new[ii-1][0] += F[ii][0]    
                    
                elif F[ii-1]>0 and abs(F[ii][0]+F[ii-1][0])*dt/zeta>2*a :
                    F_new[ii][0]=0
                    F_new[ii-1][0] += F[ii][0]                        
        F=F_new;     
        
    # if math.isnan(F):
    # F=F0*np.ones([1,nump])
    R = R + F*dt/zeta # update R & V 

    Rs[i,:,:]=R # accumulate particle positions at each step in an array Rs
    #Ws[i,:,:]=W # accumulate random forces at each step in an array Ws
    Fs[i,:,:]=F # accumulate all external forces at each step in array Fs
    timeMat[i]=i*dt # store time in each step in an array time
    if(i%100000==0):
        print(f'simulation is {i/nums*100} percent complete.') 


#%% plot
for i in range(nump):
    plt.plot(dt*np.linspace(0,nums-1,nums),Rs[:,i,:]%periodR)
    plt.xlabel('time(s)')
    plt.ylabel('position(m)')
#%% plot
for i in range(nump):
    plt.plot(dt*np.linspace(0,nums-1,nums),Rs[:,i,:]/periodR)
    plt.xlabel('time(s)')
    plt.ylabel('position(m)')
    
#%%
# Plot circles at the specified points
points=Rs[900000-1,:,:]
fig, ax = plt.subplots()

for i, point in enumerate(points):
    color = plt.cm.jet(i / len(R))  # Use a different color for each circle
    circle = plt.Circle((point[0]/periodR, 0),0.002 , color=color, fill=True)
    ax.add_patch(circle)


# Set axis limits for better visualization
ax.set_xlim(0.6, 0.8)
ax.set_ylim(-0.01, 0.01)

# Set labels and title
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_title('Circles at Points on a Straight Line')

# Show the plot
plt.show()
#%%
plt.plot(Fs[1:400000,4,:]);plt.plot(Fs[1:4000000,3,:]);plt.plot(Fs[1:4000000,1,:])
