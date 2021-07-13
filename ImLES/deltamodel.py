import numpy as np
import matplotlib.pyplot as plt

from BDIM_utils_Burgers import Kernel



class DeltaModel():

    def __init__(self, eps, grid):
        self.kernel = Kernel(eps)#
        self.eps = eps

        self.kernL = np.array([self.kernel.kernel(0, y) for y in grid])
        self.kernR = np.array([self.kernel.kernel(1, y) for y in grid])
        

    def model(self, u, spacing, nu, model="Shoeybi"):
        if model=="Shoeybi":
            edged2u_L = (u[2] - 2*u[1] + u[0])/(spacing**2)
            edged2u_R = (u[-1] - 2*u[-2] + u[-3])/(spacing**2)

            model = nu*(2*self.eps*self.kernL)*edged2u_L
            model[model.shape[0]//2:] = (nu*(2*self.eps*self.kernR)*edged2u_R)[model.shape[0]//2:]

            return model