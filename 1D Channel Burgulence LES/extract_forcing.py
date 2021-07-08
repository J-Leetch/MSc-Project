# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 17:16:23 2021

@author: james
"""


import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from BDIM_utils_Burgers import Flow1D #, exactSolution
from forcing import generate_forcing

import pickle

def get_force():
        
    return generate_forcing(1,n, A)


np.random.seed(100)

n= 2047
A = 100*np.sqrt(3)

t = 10
steps=int(t*10**4)

dt = t/steps


c=0


forcing = np.empty((steps, n+42))

for i in tqdm(range(steps)):

    if c*dt>1 or i==0:
        f = get_force()
        c=0
    forcing[i][21:2026] = f
    c+=1
    
    
pickle.dump(forcing, open("BDIM_forcing.p", "wb"))