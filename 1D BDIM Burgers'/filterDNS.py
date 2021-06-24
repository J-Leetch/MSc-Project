# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:21:11 2021

@author: james
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d
from BDIM_utils_Burgers import Grid1D, Kernel

plt.figure()

name = "DNS_Re500_255.csv"

grid = Grid1D(np.linspace(-0.05,1.05,281))

epsilon= 0.01
kern = Kernel(epsilon)
        
d = 0.5 - abs(0.5 - grid.x)

done=False
kernel = []
i=1
while not done:
    k = kern.kernel(i*grid.spacing, 0)
    if k!=0:
        kernel.append(k)
    else:
        break
    i+=1
l = len(kernel)
kernel.reverse()
kernel.append(kern.kernel(0,0))
kernel+=kernel[:-1][::-1]
kernel=np.array(kernel)



data = pd.read_csv(name, header=None)

DNS = np.zeros(281)
DNS[12:267] = data.iloc[:,1]

filtered_DNS = ndimage.convolve(DNS, kernel)*grid.spacing

plt.plot(data.iloc[:,0], data.iloc[:,1], label="DNS")

plt.plot(grid.x, filtered_DNS, label="Filtered DNS")

plt.legend()
plt.show()

np.savetxt("filtDNS_Re500_255.csv", np.hstack((grid.x.reshape(-1,1), (filtered_DNS).reshape(-1,1))), delimiter=",")

