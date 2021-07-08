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

class PeriodicQUICKGrid1D:
    
    def __init__(self, x):
               
        self.x = x
        self.n = self.x.shape[0]
        self.spacing = abs(self.x[1] - self.x[0])

        self.construct_D2()


    def QUICK(self, u):
        # self.faces = np.zeros((u.shape[0], 2))
        du = np.zeros(u.shape)

        for i in range(2, self.u.shape[0] -2):
            if u[i]>=0:
                w = 6/8 * u[i-1] + 3/8 * u[i] - 1/8 * u[i-2]
                e = 6/8 * u[i] + 3/8 * u[i+1] - 1/8 * u[i-1]
                du[i] = (e-w)/self.spacing

            if u[i]<0:
                w = 6/8 * u[i] + 3/8 * u[i-1] - 1/8 * u[i+1]
                e = 6/8 * u[i+1] + 3/8 * u[i] - 1/8 * u[i+2]
                du[i] = (e-w)/self.spacing

        return du

    def construct_D2(self):
        self.D2 = np.zeros((self.n+4, self.n+4)) #add ghost cells
                
        for i in range(2, self.D2.shape[0] -2):
            self.D2[i, i-1] = 1 / self.spacing**2
            self.D2[i, i] = -2 / self.spacing**2
            self.D2[i, i+1] = 1 / self.spacing**2



       
        
class Flow1D(PeriodicGrid1D):
    
    def __init__(self, t, n, steps, nu, model=0):
        
        super().__init__(np.linspace(0,1,n))
        self.model=model
        self.steps=steps
        self.dt = t/steps
        self.nu = nu
        
        self.u = np.zeros(self.n+2)

        self.A = 100*np.sqrt(3)

        self.i = 0
        self.c=0
        self.U = np.zeros((self.steps, self.n))

        self.INC = 100
        self.ltke = 0
        self.convprinted=False
        self.forcing=np.zeros(self.n+2)

    def f(self, u): 
        # u[0] = u[-4]
        # u[1] = u[-3]
        # u[-1] = u[3]
        # u[-2] = u[2]
        # self.force(u)
        u[0] = u[-2]
        u[-1] = u[1]

        ret = self.nu * np.matmul(self.D2, u)\
                            - np.matmul(self.D, u**2) - np.matmul(self.D, self.subgrid(u)['tau'])
        ret[1:-1] += self.forcing
        return ret #- self.QUICK(u*u))  


    def step_rk4(self):
        if self.i*self.dt>1 or self.i==0:
            self.get_force()
            self.c=0

        # self.get_force()
        # self.u[1:-1] += self.dt*self.forcing

        # self.u[0] = self.u[-4]
        # self.u[1] = self.u[-3]
        # self.u[-1] = self.u[3]
        # self.u[-2] = self.u[2]

        self.u[0] = self.u[-2]
        self.u[-1] = self.u[1]

        k1 = self.dt * self.f(self.u)
        k2 = self.dt * self.f(self.u + 0.5 * k1)
        k3 = self.dt * self.f(self.u + 0.5 * k2)
        k4 = self.dt * self.f(self.u + k3)
        self.u = self.u + (k1 + 2*(k2 + k3) + k4) / 6

        # if self.c%1000==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")

        # TKE = 0.5*self.u[1:-1]**2

        # E_k = np.fft.fft(TKE)
        self.U[self.i,:] = self.u[1:-1] # E_k
        
    
        # if abs(self.INC)>10e-3:
        #     if self.i>0 and self.i%100==0:
        #         TKE_spectrum = np.mean(np.abs(np.fft.rfft(self.E[self.i-1000:self.i]*self.spacing, axis=1)), axis=0)
        #         # RMSu = np.sqrt(np.mean(self.u**2))
        #         self.A+=self.INC
        #         print(self.A, self.INC)

        #         if (RMSu < 1 and self.RMSu_lst>1) or (RMSu > 1 and self.RMSu_lst<1):
        #             self.INC*=-0.5

        #         self.RMSu_lst = RMSu

        # else:
        #     if not self.convprinted:
        #         print(f"Converged at step {self.i}. A={self.A}")
        #         self.E = np.zeros((self.steps, self.n))
        #         self.i=0
        #         self.convprinted=True


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
        self.c+=1
    

    def step(self):
        if self.c*self.dt>1 or self.i==0:
            self.get_force()
            self.c=0

        self.u = self.u + self.dt * (self.nu * np.matmul(self.D2, self.u)\
                                            - np.matmul(self.D, self.u**2) - np.matmul(self.D, self.subgrid(self.u)['tau']) + self.forcing)
                                                            
        # if self.c%1000==0:
        #     print(f"Est. CFL {max(abs(self.u))*self.dt/self.spacing}")

        # if self.i%50==0:
        #     plt.plot(self.x[:], self.u[:], color='black', linewidth=1, marker='.')
        #     plt.ylabel("u")
        #     plt.xlabel("y")
        #     # plt.xlim(-0.02,0.1)
        #     plt.show()
            
        self.U[self.i,:] = self.u[1:-1]


        self.i+=1
        self.c+=1
        
        
    def get_force(self):
        
        self.forcing[1:-1] = generate_forcing(1,self.n, self.A)[0, :] #self.forcing[0,:] 
        # return u


    def subgrid(self, u):
        
        # no model
        if self.model==0:
            sgs = {
                'tau'   :   np.zeros(self.n+2),
                'coeff' :   0
            }
            return sgs
        
        # constant coefficient Smagorinsky
        if self.model==1:
            dudx = np.matmul(self.D, u)
            CS2   = 0.5**2
            d1    = np.abs(dudx)
            d2    = dudx
            d3    = d1*d2
            tau   = -2*CS2*(self.spacing**2)*d3
            coeff = np.sqrt(CS2)

            sgs = {
            'tau'   :   tau,
            'coeff' :   coeff
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
        
    
        