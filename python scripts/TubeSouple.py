
import numpy as np

import matplotlib.pyplot as plt



# Definition des parametres
N=2000#nombre d'element
gen=1
h=0.79 # doit être <=1 (le coefficient d'homothétie)
D0=0.02 # diamètre du tube en mètre (m)
Dmax=D0*(h**gen) #diametre maximal du tube en m

Rho=1.292 #masse volumique en kg/mcube

Mu=18.5e-6 #viscosite en kg/(m.s)

Phi=0.5 # Debit
L=0.18 #longueur du tube en m
Ax= L/(N-1) # Delta x

Alpha=(8*Rho/(np.pi)**2)*1/(Dmax**4)*(Phi**2)

Beta=(64/(np.pi))*Phi*Mu*(1/Dmax**4)*Ax

Pent=12# Pression d'entree
Ps=0#pression d'entree en Pascal
Pe=1#Pression moyenne en Pascal 
Pext=0 
 
  
def newton_raphson(f, x_guess=None, max_num_iter= 1, tolerance=1e-8, alpha=1.0, print_info=True ):
     
     
    if x_guess is None:
        x_guess = [Pent + i*(Ps-Pent)/N for i in range(N)]
    x = x_guess
 
     
    fx = f(x)
 
    
    err = np.linalg.norm(fx)
    
    
    iter = 0
 
    
     
    while err > tolerance and iter <= max_num_iter:
 
        
        x = x - alpha*np.linalg.solve(f.getJacobian(x),fx)
 
         
        fx = f(x)
 
         
        err = np.linalg.norm(fx)
 
         
        iter = iter + 1
 
         
        if print_info:
            
            print("En {0} iteration, on obtient {1}  Avec une erreur de {2}".format(iter,x,err))
            Y=np.linspace(0,L,N)
            
            plt.plot(Y,x,"bo")
            plt.ylabel('Pression (en Pa)')
            plt.xlabel('Longueur du tuple (en m)')
            plt.title('Variation de P en fonction de z avec Pref=20cmH2O')
            plt.show()
 
    " renvoie l’estimation racine, l’erreur ( norme 2 de f) de l’estimation et le nombre d’itérations que nous avons terminés"
    
    return (x, err, iter)
   
 
    
    return (x, err, iter)
class TupleSouple:
    def __init__(self):
        self.name = "TupleSouple"
 
    def numDims(self):
        return N
 
    def getJacobian(self,x):
        
        return    np.diag([1-4*((8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2)+Ax*(64/(np.pi))*x[N-1]*Mu*(1/Dmax**4))*(np.exp(-x[i]/Pe)/Pe)*(1+np.exp(-x[i]/Pe))**3 for i in range(N)]) + np.array([[ -1+4*((8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2)-Ax*(64/(np.pi))*x[N-1]*Mu*(1/Dmax**4))*np.exp(-x[i]/Pe)/Pe*(1+np.exp(-x[i]/Pe))**3 if i>k-2 and i<k else 0 for i in range(N)] for k in range(N)]) + np.array([[(1/x[N-1])*((1+np.exp(-x[j]))**4)*((64/(np.pi))*x[N-1]*Mu*(1/Dmax**4)-2*(8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2))+((1+np.exp(-x[j+1]))**4)*((64/(np.pi))*x[N-1]*Mu*(1/Dmax**4)+2*(8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2)) if i==N-1 and j<N-2 else 0 for i in range(N)] for j in range(N)]) 
        
             
 
    def __call__(self, x):
        f = np.zeros((N,))
        for i in range(0,N-1):
            
            f[i]= x[i+1]-x[i] +((1+np.exp(-x[i]/Pe))**4)*((64/(np.pi))*x[N-1]*Mu*(1/Dmax**4)*Ax-(8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2))+((1+np.exp(-x[i+1]/Pe))**4)*((64/(np.pi))*x[N-1]*Mu*(1/Dmax**4)*Ax+(8*Rho/(np.pi)**2)*1/(Dmax**4)*(x[N-1]**2))
        return f
    
        
        
 
 
def main_test():
    f = TupleSouple()
    xn = newton_raphson(f,tolerance=1e-7)
 
 
if __name__ == "__main__":
    main_test()