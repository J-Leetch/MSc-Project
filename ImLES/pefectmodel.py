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
        if C>0:
          sign=1
        else:
          sign=-1

        CS2   = sign*C**2
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
      
      plt.plot(y, Copt, label="Optimum coeff")
      plt.plot(y, ropt, label="r value", marker='x')
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
    dij_filename = "flowdata/channel_ImLES_Dij_0.0078125.p"


    target_LES_grid_res = 256

    DNS_u = pickle.load(open(filename, "rb"))
    DNS_uu = pickle.load(open(filenameuu, "rb"))
    DNS_forcing = pickle.load(open(forcing_filename, "rb"))
    dnsgrid = pickle.load(open(filenamegrid, "rb"))
    # Dij = pickle.load(open(dij_filename, "rb"))

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
    # delta_model = np.empty((steps-steps//t-2, flow.LES_U.shape[1]))
    tau_perf = np.empty((steps-steps//t-2, flow.LES_U.shape[1]))
    # smagorinsky = np.empty((n, steps-steps//t-2, flow.LES_U.shape[1]))

    # modified_smag = np.empty((n, steps-steps//t-2, flow.LES_U.shape[1]))

    C = np.array([0.75]) #np.linspace(-0.8,0.8, n)
    
    for i in tqdm(range(0, steps-steps//t-2), desc="Calculating subgrid terms for each step"):
        stepped_LES = flow.stepLES(i)
        # for c in range(C.shape[0]):
            
        #     tau = flow.smag(flow.LES_U[i], C=C[c])
            
        #     smagorinsky[c][i] = flow.BDIM(-np.matmul(flow.D, tau))
            # modified_smag[c][i] = smagorinsky[c][i]
            # modified_smag[c][i][1:4] = - modified_smag[c][i][1:4]
            # modified_smag[c][i][-4:-1] = - modified_smag[c][i][-4:-1]

        perfect_closure[i] = (flow.LES_U[i+1] - stepped_LES)/flow.dt
        # delta_model[i] = dModel.model(flow.LES_U[i], flow.spacing, flow.nu)

      
    # for i in np.random.randint(low=0, high=80000, size=5):
    #   plt.plot(y/eps, downsampled_tau_perf[i]+1000, label=r"$T_{ij} from DNS$")
    #   plt.plot(y/eps, perfect_closure[i]+1000, label="Perfect closure")
    #   # plt.plot(y, delta_model[i], label="Delta")
    #   plt.plot(y/eps, smagorinsky[0][i]+1000, label="Smagorinsky model")
    #   # plt.plot(y, modified_smag[0][i], label="mod smag")
    #   plt.yscale("log")
    #   plt.xlabel('x')
    #   plt.ylabel("Closure term")
    #   plt.legend()
    #   plt.show()

    # identidy_additional_modelling(downsampled_tau_perf, perfect_closure, y)

    # bestC = process_stats(smagorinsky, perfect_closure)

    # pickle.dump(perfect_closure, open("flowdata/perfect.p", "wb"))
    # pickle.dump(downsampled_tau_perf, open("flowdata/tau.p", "wb"))
    # pickle.dump(perfect_closure-downsampled_tau_perf, open("flowdata/diff.p", "wb"))
    pickle.dump(flow.LES_U, open("flowdata/ImLES_u.p", "wb"))

    def timevar_gif():
    
      for t in tqdm(range(0,19999,50), desc="Saving images for GIF"):
        plt.figure()
        plt.plot(y[:64]/eps, perfect_closure[t][:64], label="Perfect closure")
        plt.legend(loc='upper right')
        plt.ylim(-400,400)
        plt.xlabel(r"$y/\epsilon$")
        plt.savefig(f"images/{t}.png")

    def timevar_x0():
      t = np.linspace(0,80000,80000)
      x0 = np.empty(t.shape)
      for ti in tqdm(range(t.shape[0]), desc="Plotting"):
        x0[ti] = perfect_closure[ti][2]
        
      plt.figure()
      plt.plot(t, x0, label="Perfect closure at x=0")
      plt.legend(loc='upper right')
      plt.xlabel(r"$t$")
      plt.savefig(f"tplot.png")

    timevar_x0()