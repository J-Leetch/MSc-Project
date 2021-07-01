import pickle
import matplotlib.pyplot as plt
import numpy as np


plots = []

for name in ["periodic_4096DNS.p", "periodic_2048DNS.p",  "channel_2048DNS.p", "channel_4096DNS.p"]:
    plots.append(pickle.load(open(name, "rb")))


plt.figure()
plt.title("1D Periodic Burgulence Kinetic Energy Spectrum")
plt.xlabel("k")
plt.ylabel(r"$E(k)$")


fbins4096 = np.fft.fftfreq(4096, 1/4096)
fbins2048 = np.fft.fftfreq(2048, 1/2048)

args4096 = np.argwhere(fbins4096>0)
args2048 = np.argwhere(fbins2048>0)

args = [args4096,args2048, args2048, args4096]

plt.yscale("log")
plt.xscale("log")    

labels= ["Periodic DNS4096", "Periodic DNS2048", "Channel DNS2048", "Channel DNS4096"]
for plot, arg, label in zip(plots, args, labels):
    plt.plot(arg, plot**2, label=label)
    
x= np.array([1,1000])
plt.plot(x, x**(-2.)/500, label=r"$k^{-2} slope$", linestyle="--", color='black')

x=np.array([2,1000])
plt.plot(x, x**(-16/3.)/2, label=r"$k^{-16/3} slope$", linestyle="--", color='grey')
# plt.loglog([1,1000], np.array([1,10**-(5*3/3)]), label=r"$k^{-5/3} slope$", linestyle="--", color='green')

plt.xlim(xmin=1)
plt.ylim(ymin=10e-13, ymax=1)
plt.legend()
plt.show()