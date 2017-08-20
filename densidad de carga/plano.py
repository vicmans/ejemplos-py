import numpy as np
import math

print "Inserte los siguientes datos"
a2=float(input('Longitud del largo del cuadrito pequeno: '))
l=float(input('lado del plano: '))
v=float(input('voltaje: '))
#a2=0.2
#l=1
#v=1
a=a2/2
N=int(l/a2)
N=N**2

print "Son ",N,"parches (cuadraditos)"


epsilon=8.85*10**(-12)
iguales=(2*a/math.pi/epsilon)*np.log(1+np.sqrt(2))
area=a2**2

#Para hacer el grillado a partir de los datos
x = np.arange(a, l, 2*a)
y = np.arange(a, l, 2*a)
xx, yy = np.meshgrid(x, y)

#print "iguales", iguales
Z=np.zeros((N,N))
xf=np.append([], xx)
yf=np.append([], yy)
for m in range(0,N):
    for n in range(0,N):
      if m==n:
        Z[m,n]=iguales;
      else:
        Z[m,n]=area/np.sqrt( (xf[m]-xf[n])**2 + (yf[m]-yf[n])**2 )
print(Z)
V=np.zeros((N, 1))
for i in range(0,N):
    V[i]=v*4*math.pi*8.854187*10**(-12);

I=np.linalg.solve(Z,V);
N=int(np.sqrt(N))
I=np.reshape(I, (N, N))
#para graficar
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig = plt.figure()
ax = fig.gca(projection='3d')

colortuple = ('y', 'b')
colors = np.empty(xx.shape, dtype=str)
for y in range(len(yy)):
    for x in range(len(xx)):
        colors[x, y] = colortuple[(x + y) % len(colortuple)]

# Plot the surface with face colors taken from the array we made.
## Mejorar el color:
#my_col = cm.jet(Z/np.amax(Z))

surf = ax.plot_surface(xx, yy, I, rstride=1, cstride=1, color='b',
        linewidth=0.1, antialiased=False)

plt.title(u'Densidad de carga de un plano')
plt.xlabel('Eje x')
plt.ylabel('Eje y')
# Mostramos en pantalla
plt.show()
	
