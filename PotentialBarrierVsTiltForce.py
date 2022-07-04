# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:48:00 2022

@author: jain_
"""


k_b = 1.38e-23;

T=300;
DeltaUMax = 10*k_b*T;
R = 1e-5;


f= np.linspace(1e-17,10e-15,1000)

thetaMinU=pi-np.arcsin(-f*R/(10*k_b*T))
thetaMaxU=np.arcsin(-f*R/(10*k_b*T))+2*pi
print(thetaMinU/2/pi,thetaMaxU/2/pi)
plt.plot(f,thetaMinU/2/pi); plt.plot(f,thetaMaxU/2/pi);

#%%
DeltaU = Utilt_theta(thetaMaxU)-Utilt_theta(thetaMinU)

plt.plot(f,DeltaU/kBT)
plt.xticks(fontsize=14)
plt.xlabel('tilt force(N)',fontsize=20)
plt.ylabel('$\Delta U/k_bT$',fontsize = 20)
plt.yticks(fontsize=20)


#%%

R=10e-6        #Radius of the outer circle
periodR = 2*pi*R;
kBT = k_b*T;
F0 = 1e-17   # Drive force

rMat= np.linspace(0.0,periodR,1000);
thetaMat = np.linspace(0.0,2*pi,1000)
#Potential Function
Ufunc_ = 10*k_b*T*np.cos(thetaMat) 
#Ufunc_ = Ufunc_- max(Ufunc_)
Uforce_ = -np.diff(Ufunc_)/np.diff(rMat)


fig, axs = plt.subplots(3)
fig.suptitle('Potential, force, and tilted potential')
axs[0].plot(rMat,Ufunc_)


rMat_ = np.linspace(0.0,periodR,999)
axs[1].plot(rMat_,Uforce_)

from scipy.interpolate import interp1d

Ufunc = interp1d(rMat,Ufunc_)
Uforce = interp1d(rMat_,Uforce_)
Utilt = interp1d(rMat,Ufunc_-F0*rMat)
Utilt_theta = interp1d(thetaMat,Ufunc_-F0*rMat)

axs[2].plot(rMat/periodR,Ufunc_-F0*rMat)
