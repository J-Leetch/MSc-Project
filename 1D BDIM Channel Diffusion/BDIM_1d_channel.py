# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from BDIM_utils import Flow1D, exactSolution

nPoints = 500
eps = 0.01
t = 1
steps=int(1e5)
Re = 10000

y = np.linspace(-0.05, 1.05, nPoints)

flow = Flow1D(y, eps, t, steps, Re)

for i in tqdm(range(steps)):
    flow.step()


# print(flow.u)
plt.plot(flow.x[:], flow.u[:], color='black', linewidth=1, label=r"O2 BDIM @ $t_{final}$")

y = np.linspace(0,1,100)
# plt.plot(y, exactSolution(y, 500, Re), label="Exact solution")
plt.xlim(-0.02,0.1)
plt.ylim(0,1)
plt.ylabel("u")
plt.xlabel("y")
plt.plot([-0.01,-0.01], [0,1.1], linestyle='--', color='orange', label='Smoothing region bounds')
plt.plot([0.01,0.01], [0,1.1], linestyle='--', color='orange')
plt.legend()
plt.title(f"1D Unsteady Channel flow - Re = {Re}")
