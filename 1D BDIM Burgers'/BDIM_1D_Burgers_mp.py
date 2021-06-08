# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from BDIM_utils_Burgers import Flow1D #, exactSolution

import multiprocessing as mp
from multiprocessing import Pool

np.random.seed(100)

nPoints = 255
eps = 0.01

t = 1000000
steps=int(1e8)
Re = 500

variance = 0.01
L = 1

delta_Tr = 10

extend = 0.02

##############################


v_scale = np.sqrt(variance*L)
nu = v_scale*L/Re

y = np.linspace(-extend, L+extend, nPoints)

def func(tup):
    t, steps, pid = tup
    flow = Flow1D(y, eps, t, steps, nu, variance, delta_Tr)
    
    if pid==0:
        for i in tqdm(range(steps)):
            flow.step()
        print("Wait for all processes to finish soon...")
    else:
        for i in range(steps):
            flow.step()
        
    flow.time_average = flow.uall/steps
    return flow.time_average

if __name__=="__main__":
    nproc = mp.cpu_count()-1
    
    arg = [(t//nproc, steps//nproc, i) for i in range(nproc)]
    
    with Pool(nproc) as p:
        res = p.map(func, arg)
        
    time_average = np.sum(np.array(res), axis=0)/nproc
    
    # print(flow.u)
    plt.plot(y[:], time_average[:]/v_scale, color='black', linewidth=1, label=r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")
    
    y = np.linspace(0,1,100)
    
    plt.ylabel(r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")
    plt.xlabel("y")
    # plt.plot([-0.01,-0.01], [0,1.1], linestyle='--', color='orange', label='Smoothing region bounds')
    # plt.plot([0.01,0.01], [0,1.1], linestyle='--', color='orange')
    plt.legend()
    plt.title(f"1D Burgers' flow with Stochastic Forcing - Re = {Re}")

    plt.savefig("temp.png")
    plt.show()

    # np.savetxt("DNS_Re500_255.csv", time_average/v_scale, delimiter=",")
