# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from BDIM_utils_Burgers import Flow1D #, exactSolution

from multiprocessing import Pool
import multiprocessing as mp

np.random.seed(100)

nPoints = 4096
eps=0.0

nProc = 12

t = 1/nProc

L=1
extend=0.0

steps=int(t/(0.5*(L+extend*2)/nPoints))*4


nu = 0.001

##############################

y = np.linspace(-extend, L+extend, nPoints)


def func(i):
    flow = Flow1D(y, eps, t, steps, nu)
            
    if i==0:
        for i in tqdm(range(steps)):
            flow.step_rk4()
    else:
        for i in range(steps):
            flow.step_rk4()

    return flow.E

    
if __name__=="__main__":

    print(f"Using {nProc} processes")
    with Pool(nProc) as p:
        
        Es = p.map(func, range(nProc))

    E = np.vstack(Es)

    TKE_spectrum = np.mean(np.abs(np.fft.fft(E[:, :], axis=1)), axis=0)

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


