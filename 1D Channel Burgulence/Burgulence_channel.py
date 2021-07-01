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

nPoints = 2048
eps=0.005
t = 1

L=1
extend=0.01

steps=int(t/(0.5*(L+extend*2)/nPoints))*4


nu = 0.001

##############################

y = np.linspace(-extend, L+extend, nPoints)

flow = Flow1D(y, eps, t, steps, nu)


for i in tqdm(range(steps)):
    flow.step_rk4()


TKE_spectrum = np.mean(np.abs(np.fft.fft(flow.E[:, :], axis=1)), axis=0)

fbins = np.fft.fftfreq(nPoints, (L+extend*2)/nPoints)
args = np.argwhere(fbins>0)

plt.figure()
plt.title("1D Channel Burgulence Kinetic Energy Spectrum")
plt.xlabel("k")
plt.ylabel(r"$E(k)$")
plt.loglog(args, TKE_spectrum[args])
plt.loglog([1,args[-1]], np.array([1,1/args[-1]**2]), label=r"$k^{-2} slope$")
plt.xlim(xmin=1)
plt.legend()
plt.show()


import pickle

pickle.dump(TKE_spectrum[args], open("channel_2048DNS.p", "wb"))

# flow.time_average = flow.uall/steps

