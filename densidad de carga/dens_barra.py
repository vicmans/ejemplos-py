import numpy as np
import math

a=float(input('radio: '))
l=float(input('longitud: '))
N=input('# de elementos: ')
v=input('voltaje: ')

delta=float(l/float(N))
print delta

iguales=2*math.log((delta/2 + math.sqrt(a**2 + (delta/2)**2))/a)
print "iguales", iguales
Z=np.zeros((N,N))
for m in range(1,N+1):
    for n in range(1,N+1):
        if m==n:
            Z[m-1,n-1]=iguales;
        elif math.fabs(m-n)<=2:
           dmas=math.fabs((n-m))*delta + delta/2
           dmenos=math.fabs((n-m))*delta - delta/2
           Z[m-1,n-1]=math.log((dmas + math.sqrt((dmas)**2 + a**2))/(dmenos + math.sqrt((dmenos)**2 + a**2)))
        else:
            dmas=math.fabs((n-m))*delta + delta/2
            dmenos=math.fabs((n-m))*delta - delta/2
            Z[m-1,n-1]=math.log(dmas/dmenos)
print(Z)
V=np.zeros((N, 1))
for i in range(0,N):
    V[i]=v*4*math.pi*8.854187*10**(-12);

I=np.linalg.solve(Z,V);
print(I)

import matplotlib.pyplot as plt
x = np.arange(0., l, delta)
plt.plot(x,I)
plt.ylabel('Densidad de carga')
plt.xlabel('Longitud barra')
plt.show()
