import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d

from forcing import generate_forcing


class PeriodicGrid1D:
    
    def __init__(self, x):
               
        self.x = x
        self.n = self.x.shape[0]
        self.spacing = abs(self.x[1] - self.x[0])
        
        self.construct_D()
        self.construct_D2()
        
        
    def construct_D(self):
        self.D = np.zeros((self.n+2, self.n+2))
                
        for i in range(1, self.D.shape[0] -1):
            self.D[i, i-1] = - 1 / (2 * self.spacing)
            self.D[i, i+1] = 1 / (2 * self.spacing)

            
    def construct_D2(self):
        self.D2 = np.zeros((self.n+2, self.n+2)) #add ghost cells
                
        for i in range(1, self.D2.shape[0] -1):
            self.D2[i, i-1] = 1 / self.spacing**2
            self.D2[i, i] = -2 / self.spacing**2
            self.D2[i, i+1] = 1 / self.spacing**2
       
        
class Flow1D(PeriodicGrid1D):
    
    def __init__(self, t, n, steps, nu):
        
        super().__init__(np.linspace(0,1,n))
        
        self.steps=steps
        self.dt = t/steps
        self.nu = nu
        
        self.u = np.zeros(self.n+2)

        self.A = 0.1

        self.i = 0

        self.E = np.zeros((self.steps, self.n), dtype=complex)


    def f(self, u): 
        self.force(u)
        return (self.nu * np.matmul(self.D2, u)\
                            - np.matmul(self.D, u**2))  


    def step_rk4(self):
        if self.i%10**4==0:
            self.forcing = generate_forcing(min(self.steps-self.i,10**5), self.n, self.A)

        k1 = self.dt * self.f(self.u)
        k2 = self.dt * self.f(self.u + 0.5 * k1)
        k3 = self.dt * self.f(self.u + 0.5 * k2)
        k4 = self.dt * self.f(self.u + k3)
        self.u = self.u + (k1 + 2*(k2 + k3) + k4) / 6


        self.u[0] = self.u[-2]
        self.u[-1] = self.u[1]

        TKE = 0.5*self.u[1:-1]**2

        E_k = np.fft.fft(TKE)
        self.E[self.i,:] = E_k

        # if self.i%100==0: # and self.i>10000 and self.i<20000:
        #     plt.figure()
        #     plt.title("Periodic Burgulence - Re=1000")
        #     plt.xlim(0,1)
        #     plt.plot(self.x[:], self.u[1:-1], color='black', linewidth=1, marker='')
        #     plt.ylabel("u")
        #     plt.xlabel("y")
        #     plt.show()
        #     # plt.savefig(f"1D Burgulence Periodic/images/{self.i//100}.png")
            
        self.i+=1
    

    # def step(self):

    #     self.u = self.u + self.dt * (self.nu * np.matmul(self.D2, self.u)\
    #                                            - np.matmul(self.D, self.u**2))                                                                       
        
    #     self.force()

    #     self.u[0] = self.u[-2]
    #     self.u[-1] = self.u[1]
            
    #     if self.i%100==0:
    #         plt.figure()
    #         plt.plot(self.x[:], self.u[1:-1], color='black', linewidth=1, marker='.')
    #         plt.ylabel("u")
    #         plt.xlabel("y")
    #         plt.show()
            
    #     self.i+=1
        
        
    def force(self, u):

        u[1:-1] += self.forcing[self.i, :] #generate_forcing(1,self.n, self.A)[0, :]
        return u
        
    
        