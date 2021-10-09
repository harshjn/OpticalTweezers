# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 06:09:12 2021

@author: Tifr Anit 2
"""

# -*- coding: utf-8 -*-
""" 
Created on Tue Mar  6 10:33:02 2018

Algorithm to Calculate U(theta) given the trajectory

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt
#import timeit
pi=np.pi;
dtheta=0.0001;
phi=np.arange(0,2*pi,dtheta);
import LoadTraj as LT

filename = '1/'+F[11]
[t_data, periodX, dt]=LT.LoadData(filename);


#%
a=2e-6 # Micrometers Size of particle
r=10e-6 #Radius of the outer circle
eta=8.9e-4 #kg m^{-1}s^{-1} Dynamic Viscosity of water
k_b=1.38e-23
T=303 #Kelvin is 25degreeCelsius
gamma=3*6*pi*eta*a
Eps=2*gamma*k_b*T
dt = t_data.dt
m=4.0*pi/3*a**3*1100 #density is 1100kg/m^3
#%
#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
V_p=188e-6; #micrometers per second
f=gamma*a*V_p
#%
import scipy.io as sio
Xs=np.ravel(t_data.Xs[1:2000])

NTraj=5

plt.plot(Xs)
#%%
Xs_ = np.mod(Xs/periodX,1)
plt.plot(Xs_)
#%%

plt.hist(Xs_,bins=100);
#%%
NBins=39; # This number is better if less probably
meanXDot=np.zeros(NBins);

XMin=np.min(Xs_)
XRange=np.arange(0,1,1/NBins)
XRange = np.concatenate((XRange,[1.0]))
for i in range(NBins):
    X_=Xs_[Xs_>=XRange[i]];
    X_=X_[X_<XRange[i+1]]
    if(len(X_)>0):
#        meanThDot[i]=(theta_[len(theta_)-1]-theta_[0])/len(theta_)
        diff_X=np.diff(X_)
        diff_X=diff_X[diff_X>-0.5/NBins]
        meanXDot[i]=(np.mean(diff_X)/dt)
        
plt.plot(meanXDot)
#%%
nn=38
F0=t_data.F0
Coef_f=t_data.F0*np.ones(nn);
xMat = np.arange(0,1,1/nn)
xMat = np.concatenate((xMat,[1.0]))

import scipy.interpolate

xDot_interp = scipy.interpolate.interp1d(xMat,meanXDot)


nn=1000
Coef_f=t_data.F0*np.ones(nn);
xMat = np.arange(0,1,1/nn)
xMat = np.concatenate((xMat,[1.0]))

xDotMat=np.multiply(-1*gamma*r,xDot_interp(xMat))
xDot_int=[np.sum(xDotMat[0:i])*1/nn for i in range(nn)]

nn=999
Coef_f=t_data.F0*np.ones(nn);
xMat = np.arange(0,1,1/nn)
xMat = np.concatenate((xMat,[1.0]))

plt.plot( xMat,np.multiply(xDot_int,6.3) + (F0)*xMat+1e-13 ) 

plt.plot(np.arange(0,1,1/100),t_data.UfuncMat/r/10)
