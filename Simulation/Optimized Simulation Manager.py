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
import sys
import threading
from functools import partial

def calculate_analytical_Fc(d, alpha, r=10e-6, T=300):
    """Calculate single-particle critical force"""
    k_b = 1.38e-23
    kBT = k_b * T
    Fc = d * kBT / (alpha * r)
    return Fc

def create_robust_interpolation(thetaMat1, thetaMat2, thetaMat3, alpha, L, d, k_b, T):
    """Same interpolation function from original code"""
    epsilon = 1e-10
    
    thetaMat1 = np.linspace(0, np.pi*(1-alpha)-epsilon, 1000)
    thetaMat2 = np.linspace(np.pi*(1-alpha)+epsilon, np.pi*(1+alpha)-epsilon, 1000)
    thetaMat3 = np.linspace(np.pi*(1+alpha)+epsilon, 2*np.pi, 1000)
    
    Ufunc_2 = d*k_b*T*np.cos((thetaMat2-np.pi+np.pi*alpha)/alpha)
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
    
    Ufunc_2 = d*k_b*T*np.cos((thetaMat2-pi+pi*alpha)/alpha)
    Ufunc_1 = max(Ufunc_2)*np.ones(np.size(thetaMat1))+1e-20
    Ufunc_3 = max(Ufunc_2)*np.ones(np.size(thetaMat3))+1e-20
    
    Uforce = create_robust_interpolation(thetaMat1, thetaMat2, thetaMat3, alpha, L, d, k_b, T)
    
    # Collision forces
    def hard_sphere_force(r_ij, a):
        # r_ij = x_i - x_j (signed scalar, 1D)
        sigma    = 2.0 * a                 # contact distance
        r_cutoff = 2**(1/6) * sigma        # WCA cutoff (~1.122 * 2a)
        epsilon  = 10 * F0 *a / 12.0           # F at contact ~ F0
        F_max    = 100.0 * F0              # safety cap, 100x typical force
        
        r         = abs(r_ij)
        direction = np.sign(r_ij)
        
        if r < r_cutoff and r > 1e-10:
            sr6       = (sigma / r)**6
            sr12      = sr6**2
            magnitude = (4 * epsilon / r) * (12 * sr12 - 6 * sr6)
            magnitude = min(magnitude, F_max)   # safety cap
            return magnitude * direction
        
        return 0.0
    
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
        
        # Progress reporting
        if i % 200000 == 0 and i > 0:
            progress = i/nums*100
            if progress % 10 == 0:  # Report every 10%
                print(f'    [{os.getpid()}] {progress:.0f}% complete')
    
    # Save trajectory data to CSV
    if save_trajectories and trajectory_file and trajectory_data:
        try:
            df_traj = pd.DataFrame(trajectory_data, columns=columns)
            df_traj.to_csv(trajectory_file, index=False)
            print(f'    [{os.getpid()}] Trajectory saved: {len(trajectory_data)} timepoints → {os.path.basename(trajectory_file)}')
        except Exception as e:
            print(f'    [{os.getpid()}] Warning: Failed to save trajectory: {e}')
    
    NetDisplacement = np.mean(R - R_initial)
    TotalTime = nums*dt
    JDC = NetDisplacement / TotalTime
    
    return JDC

def simulate_single_point_with_trajectories(task_params):
    """Simulation with trajectory saving"""
    n = task_params['n']
    Fc_over_F = task_params['Fc_over_F']
    F_over_Fc = task_params['F_over_Fc']
    d = task_params['d']
    a = task_params['a']
    alpha = task_params['alpha']
    dt = task_params['dt']
    nums = task_params['nums']
    task_id = task_params['task_id']
    output_dir = task_params['output_dir']
    save_trajectories = task_params.get('save_trajectories', True)
    save_frequency = task_params.get('save_frequency', 100)
    
    # Calculate force
    Fc = calculate_analytical_Fc(d=d, alpha=alpha)
    F = F_over_Fc * Fc
    
    print(f"[{task_id:3d}] Starting: n={n:2d}, Fc/F={Fc_over_F:5.2f}, F={F*1e15:6.1f}fN", flush=True)
    
    try:
        t_start = time.time()
        
        # Set up trajectory file path
        trajectory_file = None
        if save_trajectories:
            trajectory_filename = f'trajectory_n{n:02d}_FcF{Fc_over_F:04.1f}_{task_id:03d}.csv'
            trajectory_file = os.path.join(output_dir, trajectory_filename)
        
        # Run simulation with trajectory saving
        J = CalculateJ_DC_with_trajectories(
            F, d, n, a, alpha, dt, nums, 
            save_trajectories=save_trajectories,
            trajectory_file=trajectory_file,
            save_frequency=save_frequency
        )
        
        t_end = time.time()
        
        # Calculate derived quantities
        eta = 8.9e-4
        zeta = 3 * 6 * np.pi * eta * a
        gamma_mu = zeta * J / F
        
        result = {
            'task_id': task_id,
            'n': n,
            'Fc_over_F': Fc_over_F,
            'F_over_Fc': F_over_Fc,
            'F': F,
            'Fc': Fc,
            'J': J,
            'gamma_mu': gamma_mu,
            'zeta': zeta,
            'simulation_time': t_end - t_start,
            'alpha': alpha,
            'dt': dt,
            'nums': nums,
            'save_frequency': save_frequency if save_trajectories else None,
            'trajectory_file': os.path.basename(trajectory_file) if trajectory_file else None,
            'trajectory_points': nums // save_frequency + 1 if save_trajectories else None,
            'timestamp': time.strftime("%Y%m%d_%H%M%S")
        }
        
        # Save individual result
        filename = f'result_n{n:02d}_FcF{Fc_over_F:04.1f}_{task_id:03d}.csv'
        filepath = os.path.join(output_dir, filename)
        df_single = pd.DataFrame([result])
        df_single.to_csv(filepath, index=False)
        
        print(f"[{task_id:3d}] ✓ Completed: J={J:8.2e}, Time={t_end-t_start:5.1f}s → {filename}", flush=True)
        
        return result
        
    except Exception as e:
        print(f"[{task_id:3d}] ✗ Failed: {e}", flush=True)
        error_result = {
            'task_id': task_id,
            'n': n,
            'Fc_over_F': Fc_over_F,
            'error': str(e),
            'timestamp': time.strftime("%Y%m%d_%H%M%S")
        }
        return error_result

def progress_monitor(output_dir, total_tasks, stop_event):
    """Background thread to monitor and report progress"""
    start_time = time.time()
    
    while not stop_event.is_set():
        time.sleep(30)  # Check every 30 seconds
        
        # Count completed files
        import glob
        completed_files = glob.glob(os.path.join(output_dir, "result_n*.csv"))
        trajectory_files = glob.glob(os.path.join(output_dir, "trajectory_n*.csv"))
        completed = len(completed_files)
        
        if completed > 0:
            elapsed = time.time() - start_time
            progress = completed / total_tasks
            eta_seconds = elapsed / progress - elapsed if progress > 0 else 0
            
            print(f"\n{'='*60}")
            print(f"PROGRESS UPDATE: {completed}/{total_tasks} ({progress*100:.1f}%)")
            print(f"Results: {len(completed_files)}, Trajectories: {len(trajectory_files)}")
            print(f"Elapsed: {elapsed/60:.1f} min, ETA: {eta_seconds/60:.1f} min")
            print(f"{'='*60}\n", flush=True)

def main():
    """Main optimized simulation manager with trajectory saving"""
    if len(sys.argv) < 2:
        print("Usage: python optimized_simulation_manager.py <output_directory> [--no-trajectories] [--save-freq N]")
        print("  --no-trajectories: Skip saving particle trajectories")
        print("  --save-freq N: Save trajectory every N timesteps (default: 100)")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    # Parse command line options
    save_trajectories = True
    save_frequency = 100
    
    for i, arg in enumerate(sys.argv[2:], 2):
        if arg == '--no-trajectories':
            save_trajectories = False
        elif arg == '--save-freq' and i+1 < len(sys.argv):
            save_frequency = int(sys.argv[i+1])
    
    # Read parameters from environment
    try:
        alpha = float(os.environ.get('SIMULATION_ALPHA', 0.35))
        dt = float(os.environ.get('SIMULATION_DT', 1e-2))
        d = float(os.environ.get('SIMULATION_D', 1000))
        a = float(os.environ.get('SIMULATION_A', 2e-6))
        nums = int(float(os.environ.get('SIMULATION_NUMS', 1e5)))
        n_min = int(os.environ.get('N_MIN', 1))
        n_max = int(os.environ.get('N_MAX', 10))
        max_parallel = int(os.environ.get('MAX_PARALLEL_F', 19))
    except Exception as e:
        print(f"Error reading environment variables: {e}")
        sys.exit(1)
    
    print("="*70)
    print("OPTIMIZED SIMULATION MANAGER WITH TRAJECTORY SAVING")
    print("="*70)
    print(f"Parameters:")
    print(f"  alpha={alpha}, dt={dt}, d={d}, a={a}, nums={nums}")
    print(f"  n range: {n_min} to {n_max}")
    print(f"  Max parallel F/Fc processes: {max_parallel}")
    print(f"  Save trajectories: {save_trajectories}")
    if save_trajectories:
        print(f"  Trajectory save frequency: every {save_frequency} timesteps")
        print(f"  Trajectory points per simulation: ~{nums//save_frequency + 1}")
    
    # Calculate Fc and create Fc/F range
    Fc = calculate_analytical_Fc(d=d, alpha=alpha)
    print(f"  Analytical Fc = {Fc*1e15:.1f} fN = {Fc*1e12:.3f} pN")
    
    # Fc/F values - 19 points from 0.1 to 10
    Fc_over_F_values = np.arange(1.0, 15.0+0.05, 0.1)
    F_over_Fc_values = 1.0 / Fc_over_F_values
    n_values = list(range(n_min, n_max + 1))
    
    print(f"  Fc/F values: {len(Fc_over_F_values)} points from {Fc_over_F_values.min():.1f} to {Fc_over_F_values.max():.1f}")
    print(f"  Total simulations: {len(n_values) * len(Fc_over_F_values)}")
    
    # Create all task combinations
    all_tasks = []
    task_id = 1
    
    for n in n_values:
        for i, Fc_over_F in enumerate(Fc_over_F_values):
            F_over_Fc = F_over_Fc_values[i]
            task = {
                'task_id': task_id,
                'n': n,
                'Fc_over_F': Fc_over_F,
                'F_over_Fc': F_over_Fc,
                'd': d,
                'a': a,
                'alpha': alpha,
                'dt': dt,
                'nums': nums,
                'output_dir': output_dir,
                'save_trajectories': save_trajectories,
                'save_frequency': save_frequency
            }
            all_tasks.append(task)
            task_id += 1
    
    total_tasks = len(all_tasks)
    print(f"  Total tasks created: {total_tasks}")
    
    # Estimate storage requirements
    if save_trajectories:
        trajectory_points = nums // save_frequency + 1
        avg_particles = (n_min + n_max) / 2
        bytes_per_trajectory = trajectory_points * (avg_particles + 2) * 20  # ~20 bytes per float in CSV
        total_storage_gb = total_tasks * bytes_per_trajectory / 1e9
        print(f"  Estimated trajectory storage: {total_storage_gb:.1f} GB")
    
    print("="*70)
    
    # Start progress monitoring thread
    stop_event = threading.Event()
    progress_thread = threading.Thread(
        target=progress_monitor, 
        args=(output_dir, total_tasks, stop_event)
    )
    progress_thread.daemon = True
    progress_thread.start()
    
    # Run simulations with optimal parallelization
    print(f"Starting {max_parallel} parallel processes...")
    
    t_start = time.time()
    
    with mp.Pool(processes=max_parallel) as pool:
        results = pool.map(simulate_single_point_with_trajectories, all_tasks)
    
    t_end = time.time()
    
    # Stop progress monitoring
    stop_event.set()
    
    # Process and save combined results
    print("\n" + "="*70)
    print("PROCESSING COMBINED RESULTS")
    print("="*70)
    
    # Filter successful results
    successful_results = [r for r in results if 'error' not in r]
    failed_results = [r for r in results if 'error' in r]
    
    print(f"Successful simulations: {len(successful_results)}")
    print(f"Failed simulations: {len(failed_results)}")
    
    if successful_results:
        df_all = pd.DataFrame(successful_results)
        
        # Save combined results
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        combined_file = os.path.join(output_dir, f"combined_all_results_{timestamp}.csv")
        df_all.to_csv(combined_file, index=False)
        print(f"Combined results saved: {combined_file}")
        
        # Count trajectory files
        if save_trajectories:
            import glob
            trajectory_files = glob.glob(os.path.join(output_dir, "trajectory_n*.csv"))
            print(f"Trajectory files created: {len(trajectory_files)}")
            
            # Calculate actual storage used
            total_size = sum(os.path.getsize(f) for f in trajectory_files) / 1e9
            print(f"Actual trajectory storage: {total_size:.1f} GB")
        
        # Analysis
        threshold = 1e-9
        flowing = df_all[np.abs(df_all['J']) >= threshold]
        trapped = df_all[np.abs(df_all['J']) < threshold]
        
        print(f"\nStatistics:")
        print(f"  Average simulation time: {df_all['simulation_time'].mean():.1f} ± {df_all['simulation_time'].std():.1f} seconds")
        print(f"  Total computation time: {df_all['simulation_time'].sum()/3600:.1f} hours")
        print(f"  Wall clock time: {(t_end - t_start)/60:.1f} minutes")
        print(f"  Parallel efficiency: {df_all['simulation_time'].sum()/(t_end - t_start)/max_parallel*100:.1f}%")
        print(f"  Flowing states: {len(flowing)} ({len(flowing)/len(df_all)*100:.1f}%)")
        print(f"  Trapped states: {len(trapped)} ({len(trapped)/len(df_all)*100:.1f}%)")
        
        # Save metadata
        metadata = {
            'simulation_info': {
                'total_tasks': total_tasks,
                'successful': len(successful_results),
                'failed': len(failed_results),
                'wall_time_minutes': (t_end - t_start)/60,
                'parallel_processes': max_parallel,
                'parallel_efficiency': df_all['simulation_time'].sum()/(t_end - t_start)/max_parallel,
                'trajectories_saved': save_trajectories,
                'trajectory_save_frequency': save_frequency if save_trajectories else None
            },
            'parameters': {
                'alpha': alpha, 'dt': dt, 'd': d, 'a': a, 'nums': nums,
                'n_range': [n_min, n_max],
                'F_over_Fc_range': [float(F_over_Fc_values.min()), float(F_over_Fc_values.max())],
                'Fc_analytical': Fc
            },
            'timestamp': timestamp
        }
        
        metadata_file = os.path.join(output_dir, f"simulation_metadata_{timestamp}.json")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Metadata saved: {metadata_file}")
    
    if failed_results:
        df_failed = pd.DataFrame(failed_results)
        failed_file = os.path.join(output_dir, f"failed_simulations_{timestamp}.csv")
        df_failed.to_csv(failed_file, index=False)
        print(f"Failed simulations saved: {failed_file}")
    
    print("\n" + "="*70)
    print("SIMULATION COMPLETED")
    print("="*70)
    
    return 0

if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)