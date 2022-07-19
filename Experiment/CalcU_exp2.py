# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 17:41:48 2022

@author: jain_
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# plt.style.use('ggplot')
# import sympy as sp
import math
import time
import sympy as sp
import collections
th=sp.symbols('th');
#import timeit
pi=np.pi;
#%%
#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
# import scipy.io as sio
# A=sio.loadmat('ThetaMat12Feb.mat')
# Rs = np.array(A['theta']).T
# phaseR=np.array(A['thetas']).T [0]
# plt.plot(phaseR)

Rs = np.loadtxt('thetaNonStuck2.txt')
phaseR=np.loadtxt('thetasNonStuck2.txt')+pi
plt.plot(phaseR)
len(phaseR)
#%%
#%
SizeFactor = 1
nums=len(phaseR) # Simulation steps
dt = 1/254; # Delta t of the simulation (Could be variable)
nump=1       # number of particles
# F0 = SizeFactor*5e-15# Drive force
T=300         #T #Kelvin is 25degreeCelsius
eta=8.9e-4    #eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
pi=np.pi;    
a=2e-6        # Micrometers Size of particle
r=10e-6        #Radius of the outer circle
periodR = 2*pi*r;
k_b=1.38e-23
zeta=3*6*pi*eta*a  #A factor 3 is needed to account for the plate interaction
m=4.0*pi/3*a**3*1100 #density is 1100kg/m^3
kBT = k_b*T;
#%
#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
V_p=20-6; #micrometers per second
# Initial Conditions: Potential well
amp = 10*k_b*T
Eps=2*zeta*k_b*T
std  =  np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
#%%
plt.plot(dt*np.arange(1000,2000),phaseR[1000:2000])
minX=0.6*2*pi;maxX=minX+2*pi/SizeFactor;
strtX=minX-(maxX-minX);
# plt.hlines(minX,0,2000,'r')   #Minima of potential well
# plt.hlines(maxX,0,2000,'r')   #Maxima of potential well
# plt.hlines(strtX,0,2000,'r')  #approx well start
plt.xlabel('time (s)')
#%%
plt.hist(phaseR,bins=100);
#%%
NBins=50; # This number is better if less probably
meanThDot=np.zeros(NBins);
ThMin=np.min(phaseR)
ThRange=np.arange(0,2*pi,2*pi/NBins)
ThRange = np.concatenate((ThRange,[2*pi]))
diff_phaseR=np.diff(phaseR)

for i in range(NBins):
    Th_=phaseR[phaseR>=ThRange[i]];
    Th_=Th_[Th_<ThRange[i+1]]
    
    if(len(Th_)>0):
        Mat1=phaseR>ThRange[i];
        Mat2=phaseR<ThRange[i+1];
        a=np.where(Mat1&Mat2)[0]
        b=np.diff(a)==1        
        inds=np.where(~b)[0]
        a2= np.delete(a,inds)
        a2=np.delete(a2,-1)
        diff_Th=diff_phaseR[a2]
#        meanThDot[i]=(theta_[len(theta_)-1]-theta_[0])/len(theta_)
        # diff_Th=diff_Th[diff_Th>-0.5/NBins]
        meanThDot[i]=(np.mean(diff_Th)/dt)      
        print(i, len(diff_Th))

plt.plot(meanThDot)
plt.hlines(0,0,NBins)
#%%
nn=NBins-1
# F0=t_data.F0
# Coef_f=t_data.F0*np.ones(nn);
ThMat = np.arange(0,2*pi,2*pi/nn)
ThMat = np.concatenate((ThMat,[2*pi]))
print(len(ThMat),len(meanThDot))
import scipy.interpolate
ThDot_interp = scipy.interpolate.interp1d(ThMat,meanThDot)

#%%
nn=NBins
# F0=t_data.F0
# Coef_f=t_data.F0*np.ones(nn);
ThMat = np.arange(0,2*pi,2*pi/nn)
ThMat = np.concatenate((ThMat,[2*pi]))
ThDotMat=np.multiply(-1*zeta*r,ThDot_interp(ThMat))
ThDot_int=[np.sum(ThDotMat[0:i])*1/nn for i in range(nn)]

#%% calculate f back


#%%
ThMat = np.arange(0,2*pi,2*pi/nn)
F0=8500e-15  #femtoNewton
plt.plot( ThMat[:],(np.multiply(ThDot_int,6.28)  +1.00*(F0)*ThMat[:]-1000e-18 )*1*r) 

# plt.plot(np.arange(0,2*pi,2*pi/1000),t_data.UfuncMat)
