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
nums=int(1e5) # Simulation steps
dt =1     # Delta t of the simulation (Could be variable)
nump=2       # number of particles
F0 = 300e-15# Drive force
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
dim  = 1   # system dimension (x,y,z)
std  =  np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
# np.random.seed(7) # initialize random number generator with a seed=0
R = np.zeros([nump,dim]) # array for starting & current positions    
R[1][0]=periodR*0.25;

#for i in range(nump):
#    R[i,0] = 3.72e-5; #i/nump*periodR

W = np.zeros([nump,dim]) # array to store current random forcces
F = np.zeros([nump,dim]) # array to store external force
Rs = np.zeros([nums,nump,dim]) # array to store positions at all steps
#Ws = np.zeros([nums,nump,dim]) # array to store random forces at all steps
#Fs = np.zeros([nums,nump,dim]) # array to store external forces at all steps

timeMat  = np.zeros([nums]) # an array to store time at all steps
#%%
for i in range(nums): # repeat the following operations from i=0 to nums-1
    W = std*np.random.randn(nump,dim) # generate an array of random forces
    if (R.any()>periodR or R.any()<0):
        pR = periodR*1e20;
        R__=((R*1e20)%pR)
        R_=R__/1e20;
        F = Uforce(R_) +F0#*np.ones([1,nump])
    else:
        F = Uforce(R) +F0#*np.ones([1,nump])
        
        
        
    F_Collision=np.zeros([nump,dim])
    if(abs(R[1][0]%periodR-R[0][0]%periodR)<1.05*a):

        if(R[0][0]%periodR>R[1][0]%periodR):
            F_Collision[0][0]=+4*F0
            F_Collision[1][0]=-4*F0
        else:
            F_Collision[0][0]=-4*F0
            F_Collision[1][0]=+4*F0
                
        
    # if math.isnan(F):
    #     F=F0*np.ones([1,nump])
    R = R + F*dt/zeta +W/zeta + F_Collision*dt/zeta # update R & V 

    Rs[i,:,:]=R # accumulate particle positions at each step in an array Rs
    #Ws[i,:,:]=W # accumulate random forces at each step in an array Ws
    #Fs[i,:,:]=F # accumulate all external forces at each step in array Fs
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
