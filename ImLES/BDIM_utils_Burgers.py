import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d

from forcing import generate_forcing


class Kernel:
    
    def __init__(self, epsilon: float):
        self.epsilon = epsilon
    
    def kernel(self, x: float, y: float):
        if abs(x-y) < self.epsilon:
            return (1 + np.cos(abs(x - y) * np.pi / self.epsilon))\
                    / (2*self.epsilon)
        
        else:
            return 0.0
        
    def mu0(self, d: float):
        if d <= -self.epsilon:
            return 0.0
        
        elif abs(d) < self.epsilon:
            return 0.5 * (1 + d / self.epsilon +\
                          np.sin(d * np.pi / self.epsilon) / np.pi)
            
        else:
            return 1.0
        
        
    def mu1(self, d: float):
        if abs(d) < self.epsilon:
            return self.epsilon * (0.25 - (d / (2 * self.epsilon))**2\
                - 1/(2*np.pi) * (d * np.sin(d * np.pi / self.epsilon)\
                / self.epsilon + (1 + np.cos(d * np.pi / self.epsilon))/np.pi))
        
        else:
            return 0.0
        

class Grid1D:
    
    def __init__(self, x):
               
        self.x = x
        self.n = self.x.shape[0]
        self.spacing = abs(self.x[1] - self.x[0])
        # print(self.spacing)
        
        self.construct_D()
        self.construct_D2()
        
        
    def construct_D(self):
        self.D = np.zeros((self.n, self.n))
                
        for i in range(1, self.D.shape[0] -1):
            self.D[i, i-1] = - 1 / (2 * self.spacing)
            self.D[i, i+1] = 1 / (2 * self.spacing)

            
    def construct_D2(self):
        self.D2 = np.zeros((self.n, self.n))
                
        for i in range(1, self.D2.shape[0] -1):
            self.D2[i, i-1] = 1 / self.spacing**2
            self.D2[i, i] = -2 / self.spacing**2
            self.D2[i, i+1] = 1 / self.spacing**2
       
        
class Flow1D(Grid1D):
    
    def __init__(self, y, epsilon, t, steps, nu, model=0, saveforcing=False):
        
        super().__init__(y)

        self.model=model
        self.saveforcing = saveforcing
        self.steps=steps
        self.dt = t/steps
        self.nu = nu
        self.epsilon = epsilon

        kern = Kernel(epsilon)
        
        d = 0.5 - abs(0.5 - self.x)
        
        done=False
        self.kernel = []
        i=1
        while not done:
            k = kern.kernel(i*self.spacing, 0)
            if k!=0:
                self.kernel.append(k)
            else:
                break
            i+=1
        l = len(self.kernel)
        self.kernel.reverse()
        self.kernel.append(kern.kernel(0,0))
        self.kernel+=self.kernel[:-1][::-1]
        self.kernel=np.array(self.kernel)
        
        self.mu0 = np.array([kern.mu0(di) for di in d])
        
        self.mu1 = np.array([kern.mu1(di) for di in d])

        self.interior = np.argwhere(0.5 - abs(0.5 - self.x)>0)[:,0]
        
        self.u = np.zeros(self.n)

        self.A = 100*np.sqrt(3)

        self.i = 0
        self.c=0

        self.forcing=np.zeros(self.n)
        if self.saveforcing:
            self.forcingsave = np.zeros((self.steps, self.n))
        else:
            self.forcingsave=None

        self.U = np.zeros((self.steps, self.n))


    def BDIM1(self, arg):
        return self.mu0 * arg
    
    
    def BDIM(self, arg):
        wall_n_deriv = np.matmul(self.D, arg)
        wall_n_deriv[wall_n_deriv.shape[0]//2:] = -wall_n_deriv[wall_n_deriv.shape[0]//2:]  
        #make sure the derviative has the correct direction on both walls
        return self.mu0 * arg + self.mu1 * wall_n_deriv

    def f(self, u):
        ret = self.nu * np.matmul(self.D2, u)\
                            - np.matmul(self.D, u**2) - np.matmul(self.D, self.subgrid(u)['tau'])
        ret[self.interior] += self.forcing
        return ret


    def step_rk4(self):
        if self.c*self.dt>1 or self.i==0:
            self.get_force()
            self.c=0

        # self.get_force()
        k1 = self.dt * self.f(self.u)
        k2 = self.dt * self.f(self.u + 0.5 * k1)
        k3 = self.dt * self.f(self.u + 0.5 * k2)
        k4 = self.dt * self.f(self.u + k3)
        self.u = self.BDIM(self.u + (k1 + 2*(k2 + k3) + k4) / 6)
                                                                                   
        # if self.c%100==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")
            # print(f"Est. Re {max(abs(self.u))*self.spacing/self.nu}")

        # if self.i%10==0 :# and self.i<5000:
        #     plt.figure()
        #     plt.title("Channel Burgulence BDIM")
        #     plt.plot(self.x[:], self.u[:], color='black', linewidth=1, marker='.')
        #     plt.ylabel("u")
        #     plt.xlabel("y")
        #     plt.xlim(-0.01,1.01)
        #     plt.ylim(-1.5,1.5)
        #     plt.show()

            # plt.savefig(f"images/{self.i}.png")
            
        
        self.U[self.i,:] = self.u # E_k
        

        self.i+=1
        self.c+=1


    def step(self):
        if self.c*self.dt>1 or self.i==0:
            self.get_force()
            self.c=0

        
        
        if self.saveforcing:
            self.forcingsave[self.i,:] = self.forcing

        self.u = self.BDIM(self.u + self.dt * (self.nu * np.matmul(self.D2, self.u)\
                                            - np.matmul(self.D, self.u**2) + self.forcing)) # - np.matmul(self.D, self.subgrid(self.u)['tau'])))
                                                            
        # if self.c%1000==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")

        # if self.i%10==0:
        #     plt.plot(self.x[:], self.u[:], color='black', linewidth=1, marker='.')
        #     plt.ylabel("u")
        #     plt.xlabel("y")
        #     # plt.xlim(-0.02,0.1)
        #     plt.show()
            
        self.U[self.i,:] = self.u 


        self.i+=1
        self.c+=1
    


    def get_force(self):
        
        self.forcing[self.interior] = generate_forcing(1,self.interior.shape[0], self.A)[0, :] #self.forcing[0,:] 
        
    
    def subgrid(self, u):
        
        # no model
        if self.model==0:
            sgs = {
                'tau'   :   np.zeros(self.n),
                'coeff' :   0
            }
            return sgs
        
        # constant coefficient Smagorinsky
        if self.model==1:
            dudx = np.matmul(self.D, u)
            CS2   = 0.755**2
            d1    = np.abs(dudx)
            d2    = dudx
            d3    = d1*d2
            tau   = -2*CS2*(self.spacing**2)*d3
            self.coeff = np.sqrt(CS2)

            sgs = {
            'tau'   :   tau,
            'coeff' :   self.coeff
            }
            return sgs
        
        # dynamic Smagorinsky
        # if self.model==2:
        #     uf    = utils.filterBox(u,2)
        #     uuf   = utils.filterBox(u**2,2)
        #     L11   = uuf - uf*uf
        #     dudxf = utils.filterBox(dudx,2)
        #     T     = np.abs(dudx)*dudx
        #     Tf    = utils.filterBox(T,2)   
        #     M11   = -2*(dx**2)*(4*np.abs(dudxf)*dudxf - Tf )
        #     if np.mean(M11*M11) == 0:
        #         CS2 = 0
        #     else:
        #         CS2 = np.mean(L11*M11)/np.mean(M11*M11)
        #     if CS2 < 0: 
        #         CS2 = 0
        #     d1    = utils.dealias1(np.abs(dudx),n)
        #     d2    = utils.dealias1(dudx,n)
        #     d3    = utils.dealias2(d1*d2,n)
        #     tau   = -2*CS2*(dx**2)*d3
        #     coeff = np.sqrt(CS2)
            
        #     sgs = {
        #         'tau'   :   tau,
        #         'coeff' :   coeff
        #     }
        #     return sgs
        

        # exception when none selected
        else:
            raise Exception("Please choose an SGS model in namelist.\n\
            0=no model\n\
            1=constant-coefficient Smagorinsky\n\
            2=dynamic Smagorinsky\n\
            3=dynamic Wong-Lilly\n\
            4=Deardorff 1.5-order TKE")
        
        