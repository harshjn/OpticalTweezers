#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:59:29 2026

@author: harshjn
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimized simulation manager for 32-core system with trajectory saving
Processes 19 F/Fc values in parallel, cycling through n values efficiently
Now saves particle trajectories as separate CSV files
"""

import numpy as np
import multiprocessing as mp
import pandas as pd
import time
import json
import os
from functools import partial
pi = np.pi;
k_b = 1.38e-23
T=200;

def calculate_analytical_Fc(d, alpha, r=10e-6, T=300):
    """Calculate single-particle critical force"""
    k_b = 1.38e-23
    kBT = k_b * T
    Fc = d * kBT / (2* alpha * r) # or Fc = pi*d*kbT/(alpha*L)
    return Fc

def create_robust_interpolation(thetaMat1, thetaMat2, thetaMat3, alpha, L, d, k_b, T):
    """Same interpolation function from original code"""
    epsilon = 1e-10
    
    thetaMat1 = np.linspace(0, np.pi*(1-alpha)-epsilon, 1000)
    thetaMat2 = np.linspace(np.pi*(1-alpha)+epsilon, np.pi*(1+alpha)-epsilon, 1000)
    thetaMat3 = np.linspace(np.pi*(1+alpha)+epsilon, 2*np.pi, 1000)
    
    Ufunc_2 = d/2*k_b*T*np.cos((thetaMat2-np.pi+np.pi*alpha)/alpha)
    Ufunc_1 = max(Ufunc_2)*np.ones(np.size(thetaMat1))
    Ufunc_3 = max(Ufunc_2)*np.ones(np.size(thetaMat3))
    
    Ufunc_ = np.concatenate((Ufunc_1, Ufunc_2, Ufunc_3))
    Ufunc_ = Ufunc_ - max(Ufunc_)
    thetaMat = np.concatenate((thetaMat1, thetaMat2, thetaMat3))
    rMat = thetaMat*L/(2*np.pi)
    
    dr = np.diff(rMat)
    dU = np.diff(Ufunc_)
    Uforce_ = -dU/dr
    rMat_ = rMat[:-1] + dr/2
    
    def robust_force(r_val):
        if np.isscalar(r_val):
            if r_val < rMat_[0]:
                return Uforce_[0]
            elif r_val > rMat_[-1]:
                return Uforce_[-1]
            else:
                idx = np.searchsorted(rMat_, r_val)
                r0, r1 = rMat_[idx-1:idx+1]
                f0, f1 = Uforce_[idx-1:idx+1]
                return f0 + (r_val - r0)*(f1 - f0)/(r1 - r0)
        else:
            return np.array([robust_force(r) for r in r_val])
    
    return robust_force
#%%
def CalculateJ_DC_with_trajectories(F0, d, N, a, alpha, dt, nums, save_trajectories=True, trajectory_file=None, save_frequency=100):
    """
    Optimized version with trajectory saving capability
    
    Parameters:
    -----------
    save_trajectories : bool
        Whether to save particle trajectories
    trajectory_file : str
        Path to save trajectory CSV file
    save_frequency : int
        Save every nth timestep to reduce file size (default: every 100 steps)
    """
    
    T = 300
    eta = 8.9e-4
    pi = np.pi
    r = 10e-6
    periodR = 2*pi*r
    k_b = 1.38e-23
    zeta = 3*6*pi*eta*a
    kBT = k_b*T
    L = 2*pi*r

    # Potential setup
    thetaMat1 = np.linspace(0, pi*(1-alpha), 1000)
    thetaMat2 = np.linspace(pi*(1-alpha), pi*(1+alpha), 1000)
    thetaMat3 = np.linspace(pi*(1+alpha), 2*pi, 1000)
    
    Ufunc_2 = d/2*k_b*T*np.cos((thetaMat2-pi+pi*alpha)/alpha)
    Ufunc_1 = max(Ufunc_2)*np.ones(np.size(thetaMat1))+1e-20
    Ufunc_3 = max(Ufunc_2)*np.ones(np.size(thetaMat3))+1e-20
    
    Uforce = create_robust_interpolation(thetaMat1, thetaMat2, thetaMat3, alpha, L, d, k_b, T)
    
    # Collision forces
    def hard_sphere_force(r_ij, a):
        cutoff = 2.10*a
        k_rep = F0*2
        if abs(r_ij) < cutoff:
            direction = np.sign(r_ij)
            magnitude = k_rep
            return magnitude * direction
        return 0
    
    def hard_sphere_force_multi(positions, a, periodR):
        nump = len(positions)
        F_Collision = np.zeros_like(positions)
        for i in range(nump):
            for j in range(i+1, nump):
                r_ij = (positions[j][0]) - (positions[i][0])
                if abs(r_ij) > periodR/2:
                    r_ij = periodR - abs(r_ij)
                f_hs = hard_sphere_force(r_ij, a)
                if positions[i][0]%periodR > positions[j][0]%periodR:
                    F_Collision[i][0] += f_hs
                    F_Collision[j][0] -= f_hs
                else:
                    F_Collision[i][0] -= f_hs
                    F_Collision[j][0] += f_hs
        return F_Collision
    
    # Initialize
    std = np.sqrt(2*kBT*zeta*dt)
    R = np.zeros([N, 1])
    for ii in range(N):
        R[ii][0] = L*ii/2/N
    
    R_initial = R.copy()
    
    # Initialize trajectory storage if needed
    if save_trajectories and trajectory_file:
        # Calculate number of saved timesteps
        saved_steps = nums // save_frequency + 1
        trajectory_data = []
        
        # Create header for CSV
        columns = ['timestep', 'time'] + [f'particle_{i+1}_position' for i in range(N)]
        
        # Save initial positions
        initial_row = [0, 0.0] + [R[i][0] for i in range(N)]
        trajectory_data.append(initial_row)
    
    # Main simulation loop
    for i in range(nums):
        W = std*np.random.randn(N, 1)
        R_ = R % L
        F = Uforce(R_) + F0
        F_Collision = hard_sphere_force_multi(R, a, L)
        
        R = R + F*dt/zeta + W/zeta + F_Collision*dt/zeta
        
        # Save trajectory data at specified frequency
        if save_trajectories and trajectory_file and (i + 1) % save_frequency == 0:
            current_time = (i + 1) * dt
            row = [i + 1, current_time] + [R[j][0] for j in range(N)]
            trajectory_data.append(row)
        
        report_interval = nums/100;
        # Progress reporting
        if i % report_interval == 0 and i > 0:
            progress = i/nums*100
            if progress % 10 == 0:  # Report every 10%
                print(f'    [{progress:.0f}% complete]')
    
    # Save trajectory data to CSV
    if save_trajectories and trajectory_file and trajectory_data:
        try:
            df_traj = pd.DataFrame(trajectory_data, columns=columns)
            df_traj.to_csv(trajectory_file, index=False)
            print(f'  Trajectory saved: {len(trajectory_data)} timepoints → {os.path.basename(trajectory_file)}')
        except Exception as e:
            print(f'    Warning: Failed to save trajectory: {e}')
    
    NetDisplacement = np.mean(R - R_initial)
    TotalTime = nums*dt
    JDC = NetDisplacement / TotalTime
    
    return JDC

#%%

d = 2000;
N = 3
dt = 1e-3
a = 2e-6;
nums = int(1e5)
alpha = 0.035
r = 10e-6;
Fc = calculate_analytical_Fc(d=d, alpha=alpha)
F_over_Fc = 0.1;

F0 = F_over_Fc*Fc;
output_file = f'trajectory_N{N}_alpha{alpha}_FoverFc{F_over_Fc}.csv'


JDC = CalculateJ_DC_with_trajectories(F0, d, N, a, alpha, dt, nums,save_trajectories=True,trajectory_file=output_file)
