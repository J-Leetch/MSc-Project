import pickle
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import gradient
import scipy.stats as sp
import sys
from tqdm import tqdm
import scipy.interpolate as interp
from deltamodel import DeltaModel

sys.path.append("../1D Channel Burgulence LES")

from BDIM_utils_Burgers import Flow1D

def downsampletogrid(arr, y):
    # assert(arr.shape[1]%2!=0, "Must have DNS on odd number of grid points to use this function")
    return interp.interp1d(y, arr)

# def downsample(arr, n=1):
#     for _ in range(n):
#         arr = downsample2(arr)
#     return arr

class FlowAdapted(Flow1D):

    def __init__(self, y, epsilon, t, steps, nu, model=0):
        Flow1D.__init__(self, y, epsilon, t, steps, nu, model=0)

    #new step function
    def stepLES(self, i):
        #step downsampled DNS with no model
        return self.BDIM(self.LES_U[i] + self.dt * (self.nu * np.matmul(self.D2, self.LES_U[i])\
                                            - np.matmul(self.D, self.LES_U[i]**2) + self.LES_forcing[i]))

    
    def smag(self, u,  C=0.16):
        
        # constant coefficient Smagorinsky
        dudx = np.matmul(self.D, u)
        CS2   = C**2
        d1    = np.abs(dudx)
        d2    = dudx
        d3    = d1*d2
        tau   = -2*CS2*(self.spacing**2)*d3
        # coeff = np.sqrt(CS2)

        return tau
    
def process_stats(smagorinsky, perfect_closure):
      slope = np.empty((C.shape[0], smagorinsky.shape[1]))
    #   smagorinsky = smagorinsky.reshape((smagorinsky.shape[0], smagorinsky.shape[1]*smagorinsky.shape[2]))
    #   perfect_closure = perfect_closure.reshape((smagorinsky.shape[1:]))
      
      for j in tqdm(range(C.shape[0]), desc = "Processing stats"):
    
          for i in range(perfect_closure.shape[1]):
        
            slope[j][i] = np.linalg.lstsq(smagorinsky[j][:, i][:,np.newaxis], perfect_closure[:, i])[0]
         
      
      Copt = [C[np.argmin(abs(1-np.array(slope[:, i])))] for i in range(perfect_closure.shape[1])]
    #   slopt = [slope[i][np.argmin(abs(1-np.array(slope[:, i])))] for i in range(perfect_closure.shape[1])]

      ropt = [sp.pearsonr(smagorinsky[np.argmin(abs(1-np.array(slope[:, i]))),:,  i], perfect_closure[:, i])[0] for i in range(perfect_closure.shape[1])]
      
      plt.plot(Copt, label="Optimum coeff")
      plt.plot(ropt, label="r value")
      plt.legend()
      plt.show()
    #   print("Best coefficient: ", Copt)
    # #   print("with slope: ", slopt)
    #   print("and correlation coeff.: ", ropt)
    #   plt.plot(C, slope)
    #   plt.xlabel("Cs")
    #   plt.ylabel("Slope of linear fit")
      return Copt

def identidy_additional_modelling(tau_perf, perfect_closure, y):
    diff  = (perfect_closure - tau_perf)

    meanerr = np.mean(diff**2, axis=0)
    rmserror_spatial = np.sqrt(meanerr)

    plt.scatter(y, rmserror_spatial/max(abs(rmserror_spatial)), marker='x', label="RMS error of subgrid stress")
    plt.title("ImLES additional modelling")
    plt.xlabel("x")
    plt.ylabel("normalised RMS error(perfect closure, perfect -dTij/dx)")
    plt.xlim(min(y), 0.5)
    # plt.ylim(0,1)
    plt.plot([-min(y),-min(y)], [0,1], color='grey', linestyle='--', linewidth=1)
    plt.grid()
    plt.show()


if __name__=="__main__":

    t = 30
    steps=int(t*10**4)
    nu=0.001

    #############
    eps=0.0078125


    filename = f"flowdata/channel_ImLES_target_{eps}.p"    #filename for saved DNS u at all timesteps
    filenameuu = f"flowdata/channel_ImLES_uufilt_{eps}.p"
    filenamegrid = f"flowdata/channel_ImLES_y_{eps}.p"
    forcing_filename = "flowdata/channel_ImLES_forcing_0.0078125.p"

    target_LES_grid_res = 256

    DNS_u = pickle.load(open(filename, "rb"))
    DNS_uu = pickle.load(open(filenameuu, "rb"))
    DNS_forcing = pickle.load(open(forcing_filename, "rb"))
    dnsgrid = pickle.load(open(filenamegrid, "rb"))

    L=1
    extend=eps
    y   = np.linspace(-extend, L+extend, target_LES_grid_res)

    dModel = DeltaModel(eps, y)

    downsampler = downsampletogrid(DNS_u, dnsgrid)
    downsampled_u = downsampler(y)

    tau_perf = DNS_uu - DNS_u**2
    downsampler = downsampletogrid(tau_perf, dnsgrid)
    downsampled_tau_perf = downsampler(y)[:-1]

    downsampled_tau_perf = -np.gradient(downsampled_tau_perf, abs(y[1]-y[0]), axis=1)

    downsampler = downsampletogrid(DNS_forcing, dnsgrid)
    downsampled_forcing = downsampler(y)

    # y = np.linspace(-extend, L+extend, downsampled_u.shape[1])
    flow = FlowAdapted(y, eps, t, steps, nu, model=0)
    

    flow.LES_U = downsampled_u
    del downsampled_u
    flow.U = DNS_u
    del DNS_u
    flow.LES_forcing = downsampled_forcing
    del DNS_forcing
    del downsampled_forcing

    print("Shapes of DNS and downsampled data: ", flow.U.shape, flow.LES_U.shape)

    n= 10

    perfect_closure = np.empty((steps-steps//t-2, flow.LES_U.shape[1]))
    delta_model = np.empty((steps-steps//t-2, flow.LES_U.shape[1]))
    # tau_perf = np.empty((steps-steps//t-2, flow.LES_U.shape[1]))
    smagorinsky = np.empty((n, steps-steps//t-2, flow.LES_U.shape[1]))

    C = np.array([0.75]) #np.linspace(0.7,0.8, n)
    
    for i in tqdm(range(0, steps-steps//t-2), desc="Calculating subgrid terms for each step"):
        stepped_LES = flow.stepLES(i)
        for c in range(C.shape[0]):
            
            tau = flow.smag(flow.LES_U[i], C=C[c])
            
            smagorinsky[c][i] = flow.BDIM(-np.matmul(flow.D, tau))
        perfect_closure[i] = (flow.LES_U[i+1] - stepped_LES)/flow.dt
        delta_model[i] = dModel.model(flow.LES_U[i], flow.spacing, flow.nu)
      
    for i in np.random.randint(low=0, high=80000, size=5):
      plt.plot(y, downsampled_tau_perf[i], label="Tau")
      plt.plot(y, perfect_closure[i], label="Perfect")
      plt.plot(y, delta_model[i], label="Delta")
      plt.plot(y, smagorinsky[0][i], label="smag")
      plt.legend()
      plt.show()

    identidy_additional_modelling(downsampled_tau_perf, perfect_closure, y)

    bestC = process_stats(smagorinsky, perfect_closure)