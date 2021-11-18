import math
import numpy as np
import matplotlib.pyplot as plt
petoile=5.05
dmaxt=0.02
Pext=20.0
P15=890
P0=0.0
lt=0.018
Dt=0.02
h=0.79
u=0.0000185
flux01=(P15-P0)*(math.pi*Dt**(4)*(2*h**(3))**(14)*(1-2*h**(3)))/(128*u*lt*(1-(2*h**(3))**(15)))
L = [0]*15
for i in range (0,15):
    L[i] = i
P=[]
H=[0.018]
print (flux01)
P.append(P0)
for i in range (1, 14):
    Pi=P[i-1]+(128*lt*u*2**(-i)*flux01)/(math.pi*h**(3.0*i)*Dt**(4.0))
    P.append(Pi)
    r=lt*h**(i)
    H.append(r)
rf=lt*h**(15)
H.append(rf)
P.append(P15)
print (len(P))
print (P)
plt.plot(L,P)
plt.ylabel('Valeur Pression en Pascal')
plt.xlabel("Generation de l'arbre")
plt.title("Evolution de la pression dans l'arbre")
plt.show()