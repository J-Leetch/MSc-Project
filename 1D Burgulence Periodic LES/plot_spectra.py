import pickle
import matplotlib.pyplot as plt
import numpy as np


plots = []

for name in ["periodic_2048DNS.p","periodic_1024DNS.p","periodic_512DNS.p", "periodic_256DNS.p"]:
    plots.append(pickle.load(open(name, "rb")))


plt.figure()
plt.title("1D Periodic Burgulence Kinetic Energy Spectrum")
plt.xlabel("k")
plt.ylabel(r"$E(k)$")



fbins2048 = np.fft.rfftfreq(2048, 1/2048)
fbins1024 = np.fft.rfftfreq(1024, 1/1024)
fbins512 = np.fft.rfftfreq(512, 1/512)
fbins256 = np.fft.rfftfreq(256, 1/256)


args = [fbins2048, fbins1024, fbins512, fbins256]

  

labels= ["Periodic DNS2048", "Periodic DNS1024", "Periodic DNS512", "Periodic DNS256"]
for plot, arg, label in zip(plots, args, labels):
    plt.loglog(arg, plot, label=label)
    
plt.loglog([1,10**3], np.array([1,10**-6]), label=r"$k^{-2} slope$", linestyle="--", color='grey')

plt.xlim(xmin=1)
plt.ylim(ymin=10e-13, ymax=1)
plt.legend()
plt.show()