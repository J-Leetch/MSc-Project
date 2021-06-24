# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from Burgers_periodic_utils import Flow1D #, exactSolution

np.random.seed(100)

nPoints = 2048

t = 1
steps=int(t*10**4)

nu = 0.001

##############################

flow = Flow1D(t, nPoints, steps, nu)
    
   
for i in tqdm(range(steps)):
    flow.step_rk4()


TKE_spectrum = np.mean(flow.E, axis=0)

plt.loglog(np.linspace(1,nPoints,nPoints), TKE_spectrum)
plt.show()
# flow.time_average = flow.uall/steps




