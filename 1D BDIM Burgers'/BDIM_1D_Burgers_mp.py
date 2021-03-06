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

t = 1000000
steps=int(t*100)
Re = 500

variance = 0.01

delta_Tr = 1

extend = 0.05

nPoints = 255

eps = 2*(1+extend*2)/nPoints

##############################

v_scale = np.sqrt(variance*1)
nu = v_scale*1/Re

y = np.linspace(-extend, 1+extend, nPoints)

def func(t, steps, pid):
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
    nproc = mp.cpu_count()-2
    
    arg = [(t//nproc, steps//nproc, i) for i in range(nproc)]
    
    with Pool(nproc) as p:
        res = p.starmap(func, arg)
        
    time_average = np.sum(np.array(res), axis=0)/nproc
    
    plt.plot(y[:], time_average[:]/v_scale, color='black', linewidth=1, label=r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")
        
    plt.ylabel(r"$\frac{\langle u\rangle}{\sqrt{\sigma L}}$")
    plt.xlabel("y")
    # plt.plot([-0.01,-0.01], [0,1.1], linestyle='--', color='orange', label='Smoothing region bounds')
    # plt.plot([0.01,0.01], [0,1.1], linestyle='--', color='orange')
    plt.legend()
    plt.title(f"1D Burgers' flow with Stochastic Forcing - Re = {Re}")

    plt.savefig("1D BDIM Burgers'/temp.png")
    plt.show()

    np.savetxt("1D BDIM Burgers'/test.csv", np.hstack((y.reshape(-1,1), (time_average/v_scale).reshape(-1,1))), delimiter=",")
