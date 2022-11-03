# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 16:01:52 2022

@author: jain_
"""
Th_exp = [0.      , 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,
       0.753982, 0.879646, 1.00531 , 1.13097 , 1.25664 , 1.3823  ,
       1.50796 , 1.63363 , 1.75929 , 1.88496 , 2.01062 , 2.13628 ,
       2.26195 , 2.38761 , 2.51327 , 2.63894 , 2.7646  , 2.89027 ,
       3.01593 , 3.14159 , 3.26726 , 3.39292 , 3.51858 , 3.64425 ,
       3.76991 , 3.89557 , 4.02124 , 4.1469  , 4.27257 , 4.39823 ,
       4.52389 , 4.64956 , 4.77522 , 4.90088 , 5.02655 , 5.15221 ,
       5.27788 , 5.40354 , 5.5292  , 5.65487 , 5.78053 , 5.90619 ,
       6.03186 , 6.15752,2*pi ];

Ufunc_exp = [-1.00e-20,  2.44e-19,  1.91e-19,  4.74e-20, -1.31e-19, -7.62e-19,
       -1.42e-18, -1.68e-18, -1.53e-18, -8.40e-19, -3.88e-19,  2.76e-19,
        7.65e-20, -5.82e-19, -1.11e-18, -1.34e-18, -1.23e-18, -1.09e-18,
       -9.91e-19, -1.09e-18, -1.40e-18, -1.68e-18, -1.93e-18, -1.96e-18,
       -1.81e-18, -1.85e-18, -1.90e-18, -1.84e-18, -1.90e-18, -1.86e-18,
       -1.73e-18, -1.71e-18, -1.74e-18, -2.12e-18, -2.70e-18, -3.45e-18,
       -4.16e-18, -4.51e-18, -4.57e-18, -4.28e-18, -3.80e-18, -3.21e-18,
       -2.84e-18, -2.48e-18, -1.95e-18, -1.18e-18, -5.19e-19, -3.79e-19,
       -8.52e-20, -1.07e-20,(-1e-20+ -1.07e-20)/2];

plt.plot(Th_exp,Ufunc_exp)

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# plt.style.use('ggplot')
# import sympy as sp
import math
import time
from scipy.interpolate import interp1d
#%%
SimulationTime = 1000;   #seconds
dt =1e-2                # Delta t of the simulation (Could be variable)
nums=int(SimulationTime/dt); # Simulation steps
nump=1                  # number of particles
F0 = 600e-15              # Drive force
T=300                   # T #Kelvin is 25degreeCelsius
eta=8.9e-4              # eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
pi=np.pi;    
a=2e-6                  # Micrometers Size of particle
r=10e-6                 # Radius of the outer circle
periodR = 2*pi*r;
k_b=1.38e-23
zeta=3*6*pi*eta*a       #A factor 3 is needed to account for the plate interaction
m=4.0*pi/3*a**3*1050 #density is 1100kg/m^3
kBT = k_b*T;
#%

#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
V_p=20e-6; #micrometers per second
# Initial Conditions: Potential well

rMat= np.linspace(0.0,periodR,1000);


#%%
Ufunc = interp1d(Th_exp,Ufunc_exp)

thetaMat = np.linspace(0.0,2*pi,1000)
#Potential Function
Ufunc_ = Ufunc(thetaMat) 

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
axs[2].plot(rMat,Ufunc_-F0*rMat)
#%% Calculation
dim  = 1   # system dimension (x,y,z)
std  =  np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
# np.random.seed(7) # initialize random number generator with a seed=0
R = np.zeros([nump,dim]) # array for starting & current positions    
#for i in range(nump):
#    R[i,0] = 3.72e-5; #i/nump*periodR
W = np.zeros([nump,dim]) # array to store current random forcces
F = np.zeros([nump,dim]) # array to store external force
Rs = np.zeros([nums,nump,dim]) # array to store positions at all steps
Ws = np.zeros([nums,nump,dim]) # array to store random forces at all steps
Fs = np.zeros([nums,nump,dim]) # array to store external forces at all steps
timeMat  = np.zeros([nums])    # an array to store time at all steps
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

    # if math.isnan(F):
    #     F=F0*np.ones([1,nump])
    R = R + F*dt/zeta +W/zeta # update R & V 

    Rs[i,:,:]=R # accumulate particle positions at each step in an array Rs
    Ws[i,:,:]=W # accumulate random forces at each step in an array Ws
    Fs[i,:,:]=F # accumulate all external forces at each step in array Fs
    timeMat[i]=i*dt # store time in each step in an array time
    if(i%100000==0):
        print(f'simulation is {i/nums*100} percent complete.') 
#%% plot
for i in range(nump):
    plt.plot(dt*np.linspace(0,nums-1,nums),Rs[:,i,:]/periodR)
    plt.xlabel('time(s)')
    plt.ylabel('position(m)')
#%%
rData = Rs[:,0,0]
rDot_=np.diff(rData)
#%%
tauWindow=10000;
tau = tauWindow*dt;
jMat = (rData[tauWindow:]-rData[:-tauWindow])/tau
plt.hist(jMat,40)
#%%
[V,B]=np.histogram(jMat,bins=32)
plt.plot(B[1:],V,'o')



