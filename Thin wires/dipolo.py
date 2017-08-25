import numpy as np

lamda = 1 #longitud de onda (normalizada= siempre uno)
gap=0.01*lamda # tama{\~n}o del gap de alimentacion
eta=377 # impedancia intrinseca del vacio
L=input('longitud: ')
a=.001*lamda # radio transversal de la antena
NT=input('segmentos: ')
N=NT+1 # numero de nodos
h=L/NT # paso numerico

c1=1/(2*np.sin(2*np.pi*h))
c2=-2*np.cos(2*np.pi*h)
Z=np.zeros((N,N),complex)
for m in range(1,N+1):
  for n in range(1,N+1):
    Z[m-1,n-1]=c1*(np.exp(-1j*2*np.pi*np.sqrt(((m-n-1)*h)**2+a**2))/np.sqrt(((m-n-1)*h)**2+a**2)+c2*(np.exp(-1j*2*np.pi*np.sqrt(((m-n)*h)**2+a**2))/np.sqrt(((m-n)*h)**2+a**2))+ np.exp(-1j*2*np.pi*np.sqrt(((m-n+1)*h)**2+a**2))/np.sqrt(((m-n+1)*h)**2+a**2))

V=np.zeros((N,1),complex)
if N%2==0:
  ngap=N/2
else:
  ngap=(N+1)/2
V[ngap-1]=(-1j*2*np.pi)/(eta*gap)

# Resolviendo el sistema
I=np.linalg.solve(Z,V)

print I

# graficar

import matplotlib.pyplot as plt

plt.plot(I)
plt.ylabel('Distribucion de Corriente')
plt.xlabel('segmentos')
plt.show()
