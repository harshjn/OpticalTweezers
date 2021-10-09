# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 21:23:19 2021

@author: Tifr Anit 2
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 13:24:26 2021
@author: Tifr Anit 2
"""

import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
# from fast_histogram import histogram1d
import os
import numpy.ma as ma
# Load dataset
import collections
import time
# filename = '6e-12time0.001N1n10000000T303000.npy';
#print(f'starting to build plots for {filename}')


def LoadData(filename):
    t=np.load(filename,allow_pickle=True)
    tDict = t.item()
    MyTuple = collections.namedtuple('MyTuple', tDict)
    t_data = MyTuple(**tDict)
    dt = t_data.dt;
    periodX=t_data.periodX;
    return [t_data, periodX, dt]
#%% Plot a single Traj


for F in np.arange(1e-12,8e-12,0.5e-12):
    SimulateTraj(F0=F)

os.listdir()
#%%
filename = '1/'+F[11]
[t_data, periodX, dt] = LoadData(filename);
xData = np.ravel(t_data.Xs,order='F')
plt.plot(xData[1:1000]/periodX,'.')
#%%