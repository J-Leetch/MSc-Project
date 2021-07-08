# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 12:25:30 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import pickle

flowE = pickle.load(open("periodic_2048DNS.p", "rb"))


def downsample2(arr):
    return arr[:,::2]

def plot(flowE, nPoints):
        dx = 1/nPoints
        
        TKE_spectrum2048 = np.mean(np.abs(np.fft.rfft(flowE, axis=1)), axis=0)*dx
        
        ds = downsample2(flowE)
        ds2 = downsample2(ds)
        
        TKE_spectrum1024 = np.mean(np.abs(np.fft.rfft(ds, axis=1)), axis=0)*dx*2
        TKE_spectrum512 = np.mean(np.abs(np.fft.rfft(ds2, axis=1)), axis=0)*dx*4
        
        plt.figure()
        plt.title("1D Periodic Burgulence Kinetic Energy Spectrum")
        plt.xlabel("k")
        plt.ylabel(r"$E(k)$")
        
        plt.loglog(TKE_spectrum2048**2, label="2048")
        plt.loglog(TKE_spectrum1024**2, label="1024")
        plt.loglog(TKE_spectrum512**2, label="512")
        
        
        plt.loglog([1,10**3], np.array([1,10**-6]), label=r"$k^{-2} slope$")
        
        plt.xlim(xmin=1)
        plt.ylim(ymax=1)
        plt.legend(loc='lower left')
        plt.show()
        
        
plot(flowE, 2048)