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
    
    def __init__(self, y, epsilon, t, steps, nu):
        
        super().__init__(y)
        
        self.steps=steps
        self.dt = t/steps
        self.nu = nu

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

        self.A = 1*np.sqrt(3)

        self.i = 0
        self.c=0

        self.E = np.zeros((self.steps, self.n), dtype=complex)


    def BDIM1(self, arg):
        return self.mu0 * arg
    
    
    def BDIM(self, arg):
        wall_n_deriv = np.matmul(self.D, arg)
        wall_n_deriv[wall_n_deriv.shape[0]//2:] = -wall_n_deriv[wall_n_deriv.shape[0]//2:]  
        #make sure the derviative has the correct direction on both walls
        return self.mu0 * arg + self.mu1 * wall_n_deriv

    def f(self, u):
        self.force(u)
        return self.BDIM(self.nu * np.matmul(self.D2, u)\
                            - np.matmul(self.D, u**2))  


    def step_rk4(self):
        if self.c%10**4==0:
            self.forcing = generate_forcing(min(self.steps-self.i+1,10**4+1), self.n, self.A)
            self.c=0

        k1 = self.dt * self.f(self.u)
        k2 = self.dt * self.f(self.u + 0.5 * k1)
        k3 = self.dt * self.f(self.u + 0.5 * k2)
        k4 = self.dt * self.f(self.u + k3)
        self.u = self.BDIM(self.u) + (k1 + 2*(k2 + k3) + k4) / 6
                                                                                   
        # if self.c%100==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")
            # print(f"Est. Re {max(abs(self.u))*self.spacing/self.nu}")

        if self.i%50==0 and self.i<5000:
            plt.figure()
            plt.title("Channel Burgulence BDIM")
            plt.plot(self.x[:], self.u[:], color='black', linewidth=1, marker='.')
            plt.ylabel("u")
            plt.xlabel("y")
            plt.xlim(-0.01,1.01)
            plt.ylim(-1.5,1.5)
            # plt.show()

            plt.savefig(f"images/{self.i}.png")
            
        
        TKE = 0.1*self.u**2

        # E_k = np.fft.fft(TKE)
        self.E[self.i,:] = TKE # E_k

        self.i+=1
        self.c+=1


    def step(self):
        if self.c%10**4==0:
            self.forcing = generate_forcing(min(self.steps-self.i+1,10**4+2), self.n, self.A)
            self.c=0

        self.force()
        self.u = self.BDIM(self.u + self.dt * (self.nu * np.matmul(self.D2, self.u)\
                                               - np.matmul(self.D, self.u**2)))
                                                                                   
        # if self.c%1000==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")

        # if self.i%50==0:
        #     plt.plot(self.x[:], self.u[:], color='black', linewidth=1, marker='.')
        #     plt.ylabel("u")
        #     plt.xlabel("y")
        #     # plt.xlim(-0.02,0.1)
        #     plt.show()
            
        
        TKE = 0.5*self.u**2

        # E_k = np.fft.fft(TKE)
        self.E[self.i,:] = TKE # E_k

        self.i+=1
        self.c+=1
    


        
    def force(self, u):

        u[self.interior] += self.forcing[self.c, self.interior] #generate_forcing(1,self.n, self.A)[0, :]
        return u
        
    
        