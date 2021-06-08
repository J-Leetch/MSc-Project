# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from BDIM_utils_Burgers import Flow1D #, exactSolution

np.random.seed(100)

nPoints = 200
eps = 0.0
t = 10000
steps=int(1e6)
Re = 500

variance = 0.01
L = 1

delta_Tr = 0.1

extend = 0.0

##############################


v_scale = np.sqrt(variance*L)
nu = v_scale*L/Re

y = np.linspace(-extend, L+extend, nPoints)


flow = Flow1D(y, eps, t, steps, nu, variance, delta_Tr)
    
   
for i in tqdm(range(steps)):
    flow.step()
        
flow.time_average = flow.uall/steps

# print(flow.u)
plt.plot(flow.x[:], flow.time_average[:]/v_scale, color='black', linewidth=1, label=r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")

y = np.linspace(0,1,100)

plt.ylabel(r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")
plt.xlabel("y")
# plt.plot([-0.01,-0.01], [0,1.1], linestyle='--', color='orange', label='Smoothing region bounds')
# plt.plot([0.01,0.01], [0,1.1], linestyle='--', color='orange')
plt.legend()
plt.title(f"1D Burgers' flow with Stochastic Forcing - Re = {Re}")

plt.savefig("temp.png")
