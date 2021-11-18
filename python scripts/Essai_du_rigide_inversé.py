import math

def arbre_poumon(Pe, Ps, n):
    u = 18*(10**(-6))
    pi = math.pi
    l0 = 0.18
    D0 = 0.02
    R1= (128*u*l0)/pi*(D0**4)
    S = 1
    h= 0.79
    e=2*n-1
    L = [666] * e
    for i in range (1,n+1):
        S = S + 1/(h**(3*i))*2*i
    f = (Pe-Ps)/R1*S
    F = f
    for i in range (0, 2*n-2, 2):
        L[i] = Pe
        L[i+1] = f
        Pe=Pe-(F*R1)/(2**i)*(h**(3*i))
        f = f / (2**i)
    L[2*n-2] = Ps
    return L
print(arbre_poumon(800,0,4))