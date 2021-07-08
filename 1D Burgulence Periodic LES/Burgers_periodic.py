# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from Burgers_periodic_utils import Flow1D #, exactSolution

import pickle

np.random.seed(100)

def downsample2(arr):
    return arr[:,::2]

def plot(n1, n2, n3):
         
        DNS2048 = pickle.load(open(f"periodic_{n1}DNS.p", "rb"))
        ROM256 = pickle.load(open(f"periodic_{n2}DNS.p", "rb"))
        LES256 = pickle.load(open(f"periodic_{n3}LES_MODEL1.p", "rb"))
        
        # TKEf = 0.5*np.mean(np.abs(np.fft.rfft(flow.U**2, axis=1)), axis=0)
        TKEf_DNS = 0.5*np.mean(np.abs(np.fft.rfft(DNS2048**2, axis=1)), axis=0)
        TKEf_ROM = 0.5*np.mean(np.abs(np.fft.rfft(ROM256**2, axis=1)), axis=0)
        TKEf_LES = 0.5*np.mean(np.abs(np.fft.rfft(LES256**2, axis=1)), axis=0)
        
        
        plt.figure()
        plt.title("1D Periodic Burgulence Kinetic Energy Spectrum")
        plt.xlabel("k")
        plt.ylabel(r"$E(k)$")
        
        plt.loglog((TKEf_DNS)**2,  label="FRS", color='black')
        plt.loglog(TKEf_ROM**2, label="ROM", color='red')
        plt.loglog(TKEf_LES**2, label="Smagorinsky 0.5", color='green')
        
        plt.loglog([1,10**3], np.array([1,10**-6]), linestyle='--', label=r"$k^{-2} slope$",color='Grey')
        
        plt.xlim(xmin=1)
        # plt.ylim(ymin=10e-12)
        plt.legend(loc='lower left')
        plt.show()
        
        
def save(flow, mod, nPoints, i):
    if mod==0:    
        pickle.dump(flow.U[steps//t:i], open(f"periodic_{nPoints}DNS.p", "wb"))
    else:
        pickle.dump(flow.U[steps//t:i], open(f"periodic_{nPoints}LES_MODEL{mod}.p", "wb"))
        

    
t = 10
steps=int(t*10**4)

nPoints=2048
nPoints2 = 128
nPoints3 = 128

nu = 0.001

mod=0

##############################

flow = Flow1D(t, nPoints, steps, nu, model=mod)   
# flow2 = Flow1D(t, nPoints2, steps, nu, model=mod)
# flow3 = Flow1D(t, nPoints3, steps, nu, model=1)
    
   
for i in tqdm(range(steps)):
    flow.step()
    # flow2.step_rk4()
    # flow3.step_rk4()
    if i%(steps//t)==0 and i>steps//t:
        save(flow, mod, nPoints, i)
        # save(flow2, mod, nPoints2, i)
        # save(flow3, 1, nPoints3, i)
        
        plot(nPoints, nPoints2, nPoints3)
        
save(flow, mod, nPoints, i)
# save(flow2, mod, nPoints2, i)
# save(flow3, 1, nPoints3, i)

            



