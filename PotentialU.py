# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 00:07:29 2021

@author: Tifr Anit 2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('ggplot')
import sympy as sp

pi=np.pi;
a=2e-6 # Micrometers Size of particle
r=10e-6 #Radius of the outer circle
periodX = 2*pi*r;

#%% Initial Conditions: Potential well
# Let periodicity of R be 1
x=sp.symbols('x');
amp=1e-17
func=amp*sp.cos(6*2*sp.pi*x)
# sp.plot(func,(r,0,1))
p = sp.Piecewise(
    (amp, x <= 1/3),
    (amp, x >= 0.5),
    ( func, True )
    )
# periodR = 
Ufunc=sp.lambdify(x,p)
xMat = np.linspace(0,1,100);
plt.plot(xMat,Ufunc(xMat))
q = -1*sp.diff(p,x)
Uforce = sp.lambdify(x,q);
plt.plot(xMat,Uforce(xMat)/periodX)
plt.plot(xMat,Ufunc(xMat))

# plt.plot(xMat,Uforce(xMat)/periodX)
#%%

writeMat = [xMat,Ufunc(xMat),Uforce(xMat)/periodX]
print(writeMat)
np.savetxt('r'+str(r)+'u'+str(amp)+'.csv', np.transpose(writeMat), delimiter=",")

#%%
f = -1e-12 # The external constant drift force 

plt.plot(xMat,Ufunc(xMat)+f*xMat*periodX)

# my_data = np.genfromtxt('r'+str(amp)+'.csv', delimiter=',')
