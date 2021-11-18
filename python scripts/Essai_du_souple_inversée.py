import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

h,D0,mu,L0,rho =  0.5,0.02,0.000015,0.018,0.001
Pr= 25
Pe = 12
Ps = 0
nombregen=3


#symbols
P1 = sp.Symbol('P1')
P2 = sp.Symbol('P2')
f  = sp.Symbol('f')


#fonctions
def primi(P,i):
    return sp.ln(sp.exp(4*P/Pr)+(i/3)**8)
def tube(P,i):
    return (D0*h**i)/(1+(i/3)**2*sp.exp(-P/Pr))

#Base
x  = [P1,P2,f]
fx = [0] * 3

fx[0] = (D0**4)*(h**(4*0)*Pr*0.25)*(primi(x[0],1)-primi(Pe,0))-(32*rho/(sp.pi**2))*(x[2]/2**0)**2*sp.ln(tube(x[0],1)/tube(Pe,0)) + (128*mu*(x[2]/2**0)*L0*(h**0))/sp.pi

fx[1] = (D0**4)*(h**(4*1)*Pr*0.25)*(primi(x[1],2)-primi(x[0],1))-(32*rho/(sp.pi**2))*(x[2]/2**1)**2*sp.ln(tube(x[1],2)/tube(x[0],1)) + (128*mu*(x[2]/2**1)*L0*(h**1))/sp.pi

fx[2] = (D0**4)*(h**(4*2)*Pr*0.25)*(primi(Ps,3)-primi(x[1],2))-(32*rho/(sp.pi**2))*(x[2]/2**2)**2*sp.ln(tube(Ps,3)/tube(x[1],2)) + (128*mu*(x[2]/2**2)*L0*(h**2))/sp.pi

#Jacobien
def Jacobian(v_str, f_list):
    J = sp.zeros(nombregen,nombregen)
    for i in range(0,nombregen):
        for j in range(0, nombregen):
            J[i,j] = sp.diff(f_list[i], v_str[j])
    return J


#Remplissage
def fcal(X):
    M = sp.eye(nombregen)
    for i in range (0,nombregen):
        for j in range(0,nombregen):
            M[i,j] = sp.N(Jacobian(x,fx)[i,j].subs(P1,X[0]).subs(f,X[2]).subs(P2,X[1]))
    return M

def fcal2(X):
    M = sp.Matrix([[1],[2],[3]])
    for i in range (0,nombregen):
        M[i] = sp.N(fx[i].subs(P1,X[0]).subs(f,X[2]).subs(P2,X[1]))
    return M




#Matrix
Y = sp.Matrix([[10],[9],[0.009]])
Yb= Y


#Newton Raph Method
err = 1
tolerance = 10e-4

max_num_iter=10
iter=0
while iter < max_num_iter:
    iter  = iter + 1
    H = fcal(Yb).inv()
    Yb =Yb - (H*fcal2(Yb))
print(Yb)
A = [Ps,Yb[1],Yb[0],Pe]
B = [0,1,2,3]
plt.plot(B,A)
plt.ylabel('Pression (en cmH2O)')
plt.xlabel("Generations de l'arbre")
plt.title('Variation de P en fonction des generations')
plt.show()



