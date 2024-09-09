

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 12:27:05 2024

@author: harshjn
"""

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
from scipy.interpolate import interp1d

#%%

def CalculateJ_DC(F0,Alpha):
    nums=int(1e6) # Simulation steps
    dt =1e-1     # Delta t of the simulation (Could be variable)
    nump=20       # number of particles
    # F0 = 4.1e-15# Drive force
    T=300         #T #Kelvin is 25degreeCelsius
    eta=8.9e-4    #eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
    pi=np.pi;    
    a=2e-6        # Micrometers Size of particle
    r=10e-6        #Radius of the outer circle
    periodR = 2*pi*r;
    k_b=1.38e-23
    zeta=3*6*pi*eta*a  #A factor 3 is needed to account for the plate interaction
    # m=4.0*pi/3*a**3*1100 #density is 1100kg/m^3
    kBT = k_b*T;
    #%
    #We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
    # V_p=20-6; #micrometers per second
    # Initial Conditions: Potential well
    
    rMat= np.linspace(0.0,periodR,1000);
    thetaMat = np.linspace(0.0,2*pi,1000)
    #Potential Function
    # Alhpa = 10
    Ufunc_ = Alpha*k_b*T*np.cos(thetaMat) 
    Ufunc_ = Ufunc_- min(Ufunc_)
    Uforce_ = -np.diff(Ufunc_)/np.diff(rMat)
    
    
    # fig, axs = plt.subplots(3)
    # fig.suptitle('Potential, force, and tilted potential')
    # axs[0].plot(rMat,Ufunc_/kBT)
    
    
    rMat_ = np.linspace(0.0,periodR,999)
    # axs[1].plot(rMat_,Uforce_)
    
    
    Ufunc = interp1d(rMat,Ufunc_)
    Uforce = interp1d(rMat_,Uforce_)
    
    # axs[2].plot(rMat/periodR*2*pi,(Ufunc_-F0*rMat)/kBT)
    # We calculate deltaU and the minima position
    #%---------- Calculation
    dim  = 1   # system dimension (x,y,z)
    std  =  np.sqrt(2*kBT*zeta*dt) # calculate std for \Delta W
    # np.random.seed(7) # initialize random number generator with a seed=0
    R = np.zeros([nump,dim]) # array for starting & current positions    
    
    #for i in range(nump):
    #    R[i,0] = 3.72e-5; #i/nump*periodR
    
    W = np.zeros([nump,dim]) # array to store current random forcces
    F = np.zeros([nump,dim]) # array to store external force
    Rs = np.zeros([nums,nump,dim]) # array to store positions at all steps
    #Ws = np.zeros([nums,nump,dim]) # array to store random forces at all steps
    #Fs = np.zeros([nums,nump,dim]) # array to store external forces at all steps
    JDC_Mat = np.zeros([nump,dim]);
    timeMat  = np.zeros([nums]) # an array to store time at all steps
    
    Lambda = F0*periodR/kBT;
    print(Lambda)
    print(Lambda/Alpha)
    L = periodR;
    # k = 2*pi/L;
    #%------------------------
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
        #Ws[i,:,:]=W # accumulate random forces at each step in an array Ws
        #Fs[i,:,:]=F # accumulate all external forces at each step in array Fs
        timeMat[i]=i*dt # store time in each step in an array time
        if(i%100000==0):
            print(f'simulation is {i/nums*100} percent complete.') 
    
    for i in range(nump):
        TotalTime = nums*dt;
        JDC_Mat[i] =(Rs[-1,i,0]-Rs[0,i,0]) / TotalTime;
        
    return np.mean(JDC_Mat)

#%% Calculate
import numpy as np
from multiprocessing import Pool, cpu_count

# Constants
T = 300          # Temperature in Kelvin
pi = np.pi       
r = 10e-6        # Radius of the outer circle
periodR = 2 * pi * r
k_b = 1.38e-23
kBT = k_b * T

# Define the matrices and parameters
F0Mat = np.linspace(1e-15, 10e-15, 20)*100
LambdaMat = F0Mat * periodR / kBT 
Alpha = 1000
J1 = np.zeros([len(F0Mat), 1])


# Function to run the simulation
def run_simulation(i):
    F0 = F0Mat[i]
    result = CalculateJ_DC(F0, Alpha)
    print(f'Running simulation number {i}')
    return result

# Parallel processing
if __name__ == '__main__':
    with Pool(processes=cpu_count()) as pool:
        J1 = pool.map(run_simulation, range(len(F0Mat)))

    J1 = np.array(J1).reshape(len(F0Mat), 1)
    print("All simulations completed.")

    # Saving the output J1 as a CSV file
    output_filename = 'Alpha_1000_J1.csv'
    np.savetxt(output_filename, J1, delimiter=',', header='J1', comments='')

    print(f'Results saved to {output_filename}')


# Can you parallize this code for my M2 mac to run python on 8 parallel processors?
