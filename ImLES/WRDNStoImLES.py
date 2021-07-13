import numpy as np
import matplotlib.pyplot as plt
import pickle
from numpy.core.defchararray import array

import scipy.signal as sig
from tqdm import tqdm
import scipy.interpolate as interp

from BDIM_utils_Burgers import Kernel

class AdaptedKernel(Kernel):
    
    def __init__(self, epsilon: float):
        Kernel.__init__(self, epsilon)

    def construct_kernel(self, gridspace):
        def inc_gridspace(gridspace):
            res = 0
            while True:
                yield res
                res+=gridspace

        halfkernel = []
        for y in inc_gridspace(gridspace):
            halfkernel.append(self.kernel(0,y))
            if y> self.epsilon:
                break
        kern = np.array(halfkernel[1:][::-1] + halfkernel)
        return kern

def downsample(arr, n):
    # assert(arr.shape[1]%2!=0, "Must have DNS on odd number of grid points to use this function")
    if len(arr.shape)==2:
        return arr[:,::n]
    else:
        return arr[::n]

def extend(arr, yext, lext):
    arrext = np.zeros((arr.shape[0], yext.shape[0]))
    arrext[:, lext:-lext] = arr
    return arrext

def filter(arr):
    filtered_arr = np.empty(arr.shape)
    for i in tqdm(range(filtered_arr.shape[0]), desc="Filtering"):
        filtered_arr[i] = sig.convolve(arr[i], kern, mode="same") / sum(kern)

    return filtered_arr

def downsampletogrid(arr, y):
    # assert(arr.shape[1]%2!=0, "Must have DNS on odd number of grid points to use this function")
    return interp.interp1d(y, arr)

filename = "flowdata/channel_2048DNS.p"
forcing_filename = "flowdata/channel_2048_forcing.p"


wrDNS = pickle.load(open(filename, "rb"))
# BDIM = pickle.load(open("flowdata/channel_256DNS_BDIM.p", "rb"))

y = np.linspace(0,1,wrDNS.shape[1])

targLESgrid = 256
epsilon = 1 * 2/targLESgrid 

extension=[]
yspace=abs(y[1]-y[0])
ys = -yspace
while True:
    extension.insert(0, ys)
    if ys<-epsilon:
        break
    ys -= yspace

yext = np.hstack((np.hstack((extension, y)), (1-np.array(extension))[::-1]))
# print(yext.shape)
pickle.dump(yext, open(f"flowdata/channel_ImLES_y_{epsilon}.p", "wb"))

wrDNSext = extend(wrDNS, yext, len(extension))
del wrDNS

plt.plot(yext, wrDNSext[1000])
# plt.show()

kernel = AdaptedKernel(epsilon)
kern  =  kernel.construct_kernel(yspace)

# plt.plot(kern)
# plt.show()

filtered_u = filter(wrDNSext)
pickle.dump(filtered_u, open(f"flowdata/channel_ImLES_target_{epsilon}.p", "wb"))
plt.plot(yext, filtered_u[1000])
del filtered_u
filtered_uu = filter(wrDNSext**2)
pickle.dump(filtered_uu, open(f"flowdata/channel_ImLES_uufilt_{epsilon}.p", "wb"))
del wrDNSext
del filtered_uu

DNS_force = pickle.load(open(forcing_filename, "rb"))
DNS_force_ext = extend(DNS_force, yext, len(extension))
del DNS_force
filtered_f = filter(DNS_force_ext)
pickle.dump(filtered_f, open(f"flowdata/channel_ImLES_forcing_{epsilon}.p", "wb"))
del filtered_f
del DNS_force_ext


# plt.plot(downsample(y, 8), downsample(wrDNS, 8)[1000])
plt.show()
