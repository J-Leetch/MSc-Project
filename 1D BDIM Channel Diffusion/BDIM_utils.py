# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:01:08 2021

@author: james
"""

import numpy as np
import matplotlib.pyplot as plt

class Kernel:
    
    def __init__(self, epsilon: float):
        self.epsilon = epsilon
    
    def kernel(self, x: float, y: float) -> float:
        if abs(x-y) < self.epsilon:
            return 1 + np.cos(abs(x - y) * np.pi / self.epsilon)\
                    / (2*self.epsilon)
        
        else:
            return 0.0
        
    def mu0(self, d: float) -> float:
        if abs(d) < self.epsilon:
            return 0.5 * (1 + d / self.epsilon +\
                          np.sin(d * np.pi / self.epsilon) / np.pi)
        
        elif d <= -self.epsilon:
            return 0.0
        
        else:
            return 1.0
        
        
    def mu1(self, d: float) -> float:
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
    
    def __init__(self, y, epsilon, t, steps, Re):
        super().__init__(y)
        
        self.steps=steps
        self.dt = t/steps
        
        self.nu = 1**2 / (Re * t)
        
        kern = Kernel(epsilon)
        
        d = 0.5 - abs(0.5 - self.x)
        
        self.mu0 = np.array([kern.mu0(di) for di in d])
        
        self.mu1 = np.array([kern.mu1(di) for di in d])
                    
        self.u = np.zeros(self.n)
        for i in range(self.u.shape[0]):
            if self.x[i] > 0 and self.x[i] < 1:
                self.u[i] = 1
        
        self.i=0
        
        # plt.plot(self.x[-15:], self.mu0[-15:])
        # plt.show()
        
        
    def BDIM1(self, arg):
        return self.mu0 * arg
    
    
    def BDIM(self, arg):
        return self.mu0 * arg + self.mu1 * np.matmul(self.D, arg)
        
    
    def step(self):
        
        # if self.i%20000==0:
        #     if self.i==0:
        #         plt.plot(self.x[:], self.u[:], color='black', linewidth=1, alpha=0.2+(self.i/self.steps)*0.8, label=r"O2 BDIM @ $t_{initial}$")
        #     else:
        #         plt.plot(self.x[:], self.u[:], color='black', linewidth=1, alpha=0.2+(self.i/self.steps)*0.8)
        #     plt.show()
                    
        self.u = self.BDIM(self.u + self.dt * self.nu * np.matmul(self.D2, self.u))
        self.du = np.matmul(self.D, self.u)
                    
        self.i+=1


def exactSolution(y, k, Re):
    uExact = np.zeros(y.shape[0])
    U = 1
    
    for ki in range(k):
        uExact += np.exp(-(2*ki+1)**2 * np.pi**2 / Re) * np.sin((2*ki+1) * np.pi * y/1)
         
    
    return uExact * U
        
        
    
        