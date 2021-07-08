# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from BDIM_utils_Burgers import Flow1D #, exactSolution

import pickle

np.random.seed(100)

##############################

def plot(n1, n2, n3):
        if eps==0:
            bdim = ""
        else:
            bdim="_BDIM"
        DNS2048 = pickle.load(open(f"flowdata/channel_{n1}DNS{bdim}.p", "rb"))
        ROM256 = pickle.load(open(f"flowdata/channel_{n2}DNS{bdim}.p", "rb"))
        LES256 = pickle.load(open(f"flowdata/channel_{n3}LES{bdim}_MODEL1.p", "rb"))
        
        TKEf_DNS = 0.5*np.mean(np.abs(np.fft.rfft(DNS2048**2, axis=1)), axis=0)
        TKEf_ROM = 0.5*np.mean(np.abs(np.fft.rfft(ROM256**2, axis=1)), axis=0)
        TKEf_LES = 0.5*np.mean(np.abs(np.fft.rfft(LES256**2, axis=1)), axis=0)
        
        fbins1 = np.fft.rfftfreq(n1, 1/n1)
        fbins2 = np.fft.rfftfreq(n2, 1/n2)
        fbins3 = np.fft.rfftfreq(n3, 1/n3)
        
        
        plt.figure()
        plt.title(f"1D Channel {bdim[1:]} Burgulence Kinetic Energy Spectrum")
        plt.xlabel("k")
        plt.ylabel(r"$E(k)$")
        
        plt.loglog(fbins1, (TKEf_DNS)**2,  label="FRS", color='black')
        plt.loglog(fbins2,TKEf_ROM**2, label="ROM", color='red')
        plt.loglog(fbins3, TKEf_LES**2, label=f"Smagorinsky optimum ({flow3.coeff})", color='green')
        
        plt.loglog([1,10**3], np.array([1,10**-6])*1000, linestyle='--', label=r"$k^{-2} slope$",color='Grey')
        
        plt.xlim(xmin=1, xmax=200)
        plt.ylim(ymin=10e-4)
        plt.legend(loc='lower left')
        plt.show()
        
        
def save(flow, mod, nPoints, i):
    if eps==0:
        bdim = ""
    else:
        bdim="_BDIM"
    if mod==0:
        pickle.dump(flow.U[steps//t:i], open(f"flowdata/channel_{nPoints}DNS{bdim}.p", "wb"))  #cut out first, developing timesteps
    else:
        pickle.dump(flow.U[steps//t:i], open(f"flowdata/channel_{nPoints}LES{bdim}_MODEL{mod}.p", "wb"))
        
    try:
        pickle.dump(flow.forcingsave[steps//t:i], open(f"flowdata/channel_{nPoints}_forcing{bdim}.p", "wb"))
    except:
        pass
        

    
t = 10
steps=int(t*10**4)

nPoints= 2049
nPoints2 = 257
nPoints3 = 257

nu = 0.001

mod=0

eps=0.00

L=1
extend=0.0


y1 = np.linspace(-extend, L+extend, nPoints)
y2 = np.linspace(-extend, L+extend, nPoints2)
y3 = np.linspace(-extend, L+extend, nPoints3)



##############################

# flow = Flow1D(y1, eps, t,  steps, nu, model=mod)   
# flow2 = Flow1D(y2, eps, t,  steps, nu, model=mod, saveforcing=True)
flow3 = Flow1D(y3, eps, t,  steps, nu, model=1)
    
   
for i in tqdm(range(steps)):
    # flow.step()
    # flow2.step()
    flow3.step()
    
    # if i%(steps//t)==0 and i>steps//t:
    #     save(flow, mod, nPoints, i)
    #     # save(flow2, mod, nPoints2, i)
    #     # save(flow3, 1, nPoints3, i)
        
    #     plot(nPoints, nPoints2, nPoints3)
        
# save(flow, mod, nPoints, i)
# save(flow2, mod, nPoints2, i)
save(flow3, 1, nPoints3, i)

plot(nPoints, nPoints2, nPoints3)

