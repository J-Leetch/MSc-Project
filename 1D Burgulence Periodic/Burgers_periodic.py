# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from Burgers_periodic_utils import Flow1D #, exactSolution

# np.random.seed(100)

nPoints = 4096

t = 1
# steps=int(t*10**4)*8
steps= int(t/(0.5/nPoints))

nu = 0.001

##############################

flow = Flow1D(t, nPoints, steps, nu)

def plot(flow, nPoints, i):
    TKE_spectrum = np.mean(np.abs(np.fft.fft(flow.E[:i,128:-128], axis=1)), axis=0)*flow.spacing

    fbins = np.fft.fftfreq(nPoints, 1/nPoints)
    args = np.argwhere(fbins>0)
    
    plt.figure()
    plt.title("1D Periodic Burgulence Kinetic Energy Spectrum")
    plt.xlabel("k")
    plt.ylabel(r"$E(k)$")
    plt.loglog(args, TKE_spectrum[args])
    plt.loglog([1,args[-1]], np.array([1,1/args[-1]**2]), label=r"$k^{-2} slope$")
    plt.xlim(xmin=1)
    plt.ylim(ymin=10e-12, ymax=1)
    plt.show()
    
   
for i in tqdm(range(steps)):
    flow.step_rk4()
    if i%10000==0 and i!=0:
        plot(flow, nPoints, i)

plot(flow, nPoints, i)

# import pickle

# pickle.dump(TKE_spectrum[args], open("periodic_1024DNS.p", "wb"))

# flow.time_average = flow.uall/steps


