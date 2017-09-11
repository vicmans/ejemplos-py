""" Calculo de distribucion de corriente de un dipolo simple
con exitacion Voltaje GAP a 1V

"""

import numpy as np

lamda = 1 #longitud de onda (normalizada= siempre uno)
gap=0.01*lamda # tamano del gap de alimentacion
eta=377 # impedancia intrinseca del vacio
L=float(input('longitud (en longitudes de onda): '))
a=.001*lamda # radio transversal de la antena
NT=input('Numero de segmentos: ')
N=NT+1 # numero de nodos
h=L/N # paso numerico

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

# distribucion de corriente con funciones base
Nx=10*N # numero de muestras de la corriente
hx=L/Nx # paso para el computo de las funciones bases
x2=np.linspace(0,L,Nx) # subespacio x'
i=np.zeros((1,Nx)) # vector de muestras de la corriente
f1=np.zeros((N,Nx))
f2=np.zeros((N,Nx))

c3=1/np.sin(2*np.pi*h)

for n in range(1,N):
  uno=((x2>(n*h)) & (x2<(h*(n+1)))).astype(int)
  f1[n-1]=np.sin(2*np.pi*((n+1)*h-x2))*uno
  f2[n-1]=np.sin(2*np.pi*(x2-(n-1)*h))*(((x2<(n*h))&(x2>(h*(n-1))))).astype(int)
  i=i+I[n-1]*c3*(f1[n-1]+f2[n-1])


# vector de radiacion
y=np.arange(-L/2,L/2-hx/2,hx)#[:, newaxis]
theta=np.linspace(0,np.pi,40)
ctheta=np.diag(np.cos(theta))
xtheta=np.tile(y,(len(ctheta),1))
ytheta=xtheta.T.dot(ctheta)
Nr=i.dot(np.exp(1j*2*np.pi*ytheta)) # vector de radiacion
Atheta=Nr*np.sin(theta); # Atheta en la zona lejana
Afi=Atheta[0][19] # El valor cuando tita=90, plano xy
fi=np.linspace(0,2*np.pi,40)
# graficar

import matplotlib.pyplot as plt

fig = plt.figure()
fig.suptitle(r'Dipolo '+str(L)+'$\lambda$',fontsize=14)

ax1 = fig.add_subplot(221)
ax1.plot(np.absolute(I));
plt.ylabel('Distribucion de Corriente')
plt.xlabel('muestras')

ax2 = fig.add_subplot(223)
ax2.plot(x2,np.absolute(i[0]),color='green');
plt.ylabel('Distribucion de Corriente')
plt.xlabel('landas')

ax3 = fig.add_subplot(222, projection='polar')

F=np.absolute(Atheta[0])/np.max(np.absolute(Atheta[0]))
ax3.plot(theta.T,F,-theta.T,F, color='blue')
ax3.set_title("Patron de radiacion", va='baseline')

ax4 = fig.add_subplot(224, projection='polar')

Ffi=np.absolute(Afi/np.max(np.absolute(Atheta[0]))*np.ones((1,40)))
ax4.plot(fi,Ffi[0], color='blue')
ax4.set_title("Plano Horizontal", va='baseline')

fig.subplots_adjust(hspace=0.4)

plt.show()
