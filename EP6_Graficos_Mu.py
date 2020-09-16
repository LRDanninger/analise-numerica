
import numpy as np
from matplotlib import pyplot as plt


def functionk(x):
   sigma=0.2
   func = 2*((np.exp(1))**(((-1*(x-ComprimentoL/2)**2))/(sigma)**2))
   if abs(func)<10e-6:
   	func=0
   return func;
def functionkderivada(x):
   sigma=0.1
   func = 2*((np.exp(1))**(((-1*(x-ComprimentoL/2)**2))/(sigma)**2))*((2*x-ComprimentoL)/(sigma)**2)
   return func;
n=99
h=1/(n+1)
testeResultado = np.arange(h,1,h/10)
Plot=np.zeros(len(testeResultado))
Plot1=np.zeros(len(testeResultado))
ComprimentoL=1
for i in range (0,len(testeResultado)):
	#Plot[i]=functionkderivada(testeResultado[i])
	Plot1[i]=functionk(testeResultado[i])

print(Plot1)
#plt.plot(testeResultado,Plot)
plt.plot(testeResultado,Plot1)
plt.xlabel("Valores de X")
plt.ylabel("Função Mu")
plt.title("Função Mu")
plt.show()