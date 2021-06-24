
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from Burgers_periodic_utils import Flow1D #, exactSolution

from multiprocessing import Pool
import multiprocessing as mp

np.random.seed(100)

nPoints = 512

t = 1 /14
steps=int(t*10**4)

nu = 0.001

def func(i):
    flow = Flow1D(t, nPoints, steps, nu)
            
    if i==0:
        for i in tqdm(range(steps)):
            flow.step_rk4()
    else:
        for i in range(steps):
            flow.step_rk4()

    return flow.E
##############################


if __name__=="__main__":
    nproc = 14# mp.cpu_count()
    print(f"Using {nproc} processes")
    with Pool(nproc) as p:
        
        Es = p.map(func, range(nproc))

    E = np.vstack(Es)

    TKE_spectrum = np.mean(E, axis=0)

    plt.plot(np.linspace(1,nPoints,nPoints),TKE_spectrum)
    plt.show()
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("Energy spectrum.png")
        # flow.time_average = flow.uall/steps




