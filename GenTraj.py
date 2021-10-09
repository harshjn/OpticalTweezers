# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:54:20 2021

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('ggplot')
# import sympy as sp
import math
import collections
import time
#%% set values of physical constants of the problem
def SimulateTraj(nums=int(1e6),dt =5e-3,nump=1,F0 = 7.5e-12, T=303,eta=8.9e-4,ForceFile='r1e-05u1e-17.csv'):
    #%%
    pi=np.pi;
    a=2e-6 # Micrometers Size of particle
    r=10e-6 #Radius of the outer circle
    periodX = 2*pi*r;
    #eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
    k_b=1.38e-23
     #T #Kelvin is 25degreeCelsius
    zeta=3*6*pi*eta*a
    m=4.0*pi/3*a**3*1100 #density is 1100kg/m^3
    kBT = k_b*T;
    #%%
    #We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
    V_p=188e-6; #micrometers per second
    #%% Initial Conditions: Potential well
    # Let periodicity of R be 1
    # We load data for the U(x) potential in which our brownian particle moves
    # A = np.load(ForceFile)
    d = np.genfromtxt(ForceFile, delimiter=',')
    xMat= d[:,0];
    xMat_=xMat;
    Uforce_ = d[:,2]
    Ufunc_ = d[:,1]
    from scipy.interpolate import interp1d
    Uforce = interp1d(xMat_*periodX,Uforce_)
    Ufunc = interp1d(xMat_*periodX,Ufunc_)
    #%% Calculation
    #periodX = 1 #Set period of R as 1
    dim  = 1   # system dimension (x,y,z)
    # nump = 1 # number of independent Brownian particle trajectories
    # nums = int(1e6) # number of simulation steps
    # dt   = 0.00005 # set time increment, \Delta t
    # zeta = 1.0  # set friction constant, \zeta
    #m    = 1.0  # set particle mass, m
    #kBT  = 1.0  # set temperatute, k_B T
    # F0   = 70e-12  # set external drift force
    std  = np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
    # np.random.seed(7) # initialize random number generator with a seed=0
    X = np.zeros([nump,dim])# array for starting & current positions    
    for i in range(nump):
        X[i,0]=i/nump*periodX
    W = np.zeros([nump,dim]) # array to store current random forcces
    F = np.zeros([nump,dim]) # array to store external force
    Xs = np.zeros([nums,nump,dim]) # array to store positions at all steps
    Ws = np.zeros([nums,nump,dim]) # array to store random forces at all steps
    Fs = np.zeros([nums,nump,dim]) # array to store external forces at all steps
    
    timeMat  = np.zeros([nums]) # an array to store time at all steps
    #%%
    for i in range(nums): # repeat the following operations from i=0 to nums-1
        W = std*np.random.randn(nump,dim) # generate an array of random forces
        if (X.any()>periodX or X.any()<0):
            pX = periodX*1e20;
            X__=((X*1e20)%pX)
            X_=X__/1e20;
            F = Uforce(X_) +F0#*np.ones([1,nump])
        else:
            F = Uforce(X) +F0#*np.ones([1,nump])

        # if math.isnan(F):
        #     F=F0*np.ones([1,nump])
        X = X + F*dt/zeta +W/zeta # update R & V 
    
        Xs[i,:,:]=X # accumulate particle positions at each step in an array Xs
        Ws[i,:,:]=W # accumulate random forces at each step in an array Ws
        Fs[i,:,:]=F # accumulate all external forces at each step in array Fs
        timeMat[i]=i*dt # store time in each step in an array time
        if(i%100000==0):
            print(f'simulation is {i/nums*100} percent complete.') 
    # plt.plot(np.ravel(Xs,order='F')/periodX)            
    #%%
    Traj = collections.namedtuple('Traj', ['dt','timeMat','Xs','T','std','periodX','UforceMat','UfuncMat','xMat','eta','zeta','F0'])
    t=Traj(dt,timeMat, Xs,T,std,periodX,Uforce_,Ufunc_,xMat,eta,zeta,F0)
    np.save( str(F0)+'time'+str(dt)+'N'+str(nump)+'n'+str(nums)+'T'+str(T)+'.npy',t._asdict() ,allow_pickle =True)
    return
