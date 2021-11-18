# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:22:21 2019

@author: BIKI FOFE
"""

import numpy as np
import matplotlib.pyplot as plt

# Definition des parametres
N=1000#nombre d'element

Dmax=0.025 #diametre du tube en m

Rho=1.292e-3 #masse volumique en kg/mcube

Mu=18.5e-6 #viscosite en kg/(m.s)

Phi=0.5 # Debit
L=0.15 #longueur du tube en m
Ax= L/(N-1) # Delta x

Alpha=(8*Rho/(np.pi)**2)*1/(Dmax**4)*(Phi**2)

Beta=(64/(np.pi))*Phi*Mu*(1/Dmax**4)*Ax

Pent=1# Pression d'entree
Ps=0#pression d'entree en Pascal
Pe=12#Pression moyenne en Pascal 
Pext=0 
def newton_raphson(f, x_guess=None, max_num_iter= 1000, tolerance=1e-8, alpha=1.0, print_info=True ):
     
     
    if x_guess is None:
        x_guess = [Pent + i*(Ps-Pent)/N for i in range(N)]
    x = x_guess
 
     
    fx = f(x)
 
    
    err = np.linalg.norm(fx)
    
    
    iter = 0
 
    
     
    while err > tolerance and iter < max_num_iter:
 
        
        x = x - alpha*np.linalg.solve(f.getJacobian(x),fx)
 
         
        fx = f(x)
 
         
        err = np.linalg.norm(fx)
 
         
        iter = iter + 1
 
 
        
    if print_info:
            
            print("En {0} iteration, on obtient {1}  Avec une erreur de {2}".format(iter,x,err))
             
            Y=np.linspace(0,L,N)
            
            plt.plot(Y,x,"go")
            plt.ylabel('Pression (en Pa)')
            plt.xlabel('Longueur du tuple (en m)')
            plt.title('Variation de P en fonction de z avec Pe = 1cmH2O')
            plt.show()
 
 
    " renvoie l’estimation racine, l’erreur ( norme 2 de f) de l’estimation et le nombre d’itérations que nous avons terminés"
    
    return (x, err, iter)
   
   
class TupleRigide:
    def __init__(self):
        self.name = "TupleRigide"
 
    def numDims(self):
        return N
 
    def getJacobian(self,x):
        
        return    np.diag([-1 for i in range(N)]) + np.array([[1 if i>k and i<k+2 else 0 for i in range(N)] for k in range(N)]) 
        
             
 
    def __call__(self, x):
        f = np.zeros((N,))
        for i in range(0,N-1):
            
            f[i]= x[i+1]-x[i] + 2*Beta
        return f
    
        
        
 
 
def main_test():
    f = TupleRigide()
    xn = newton_raphson(f,tolerance=1e-7)
 
 
if __name__ == "__main__":
    main_test()
    
