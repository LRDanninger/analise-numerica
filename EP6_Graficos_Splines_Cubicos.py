
import numpy as np
from matplotlib import pyplot as plt

def functionB(x):
	if x<-2:
		func=0
	elif x>=-2 and x<=-1:
		func=(1/4)*(2+x)**3
	elif x>-1 and x<=0:
		func=(1/4)*(((2+x)**3)-4*((1+x)**3))
	elif x>0 and x<=1:
		func=(1/4)*(((2-x)**3)-4*((1-x)**3))
	elif x>1 and x<=2:
		func=(1/4)*(2-x)**3
	else:
		func=0
   	
	return func;

def functionBDerivada(x):
	if x<-2:
		func=0
	elif x>=-2 and x<=-1:
		func=(3/4)*(x+2)**2
	elif x>-1 and x<=0:
		func=(-3/4)*x*(3*x+4)
	elif x>0 and x<=1:
		func=(3/4)*x*(3*x-4)
	elif x>1 and x<=2:
		func=(-3/4)*(x-2)**2
	else:
		func=0
		
	return func;


def functionPhi(x,j):
	if j==0:
		phi=functionB((x-0*h)/h)-4*functionB((x-(-1)*h)/h)
	elif j==1:
			phi=functionB((x-1*h)/h)-functionB((x-(-1)*h)/h)
	elif j>=2 and j<=n-1:
		phi=functionB((x-j*h)/h)
	elif j==n:
		phi=functionB((x-n*h)/h)-functionB((x-(n+2)*h)/h)
	elif j==n+1:
		phi=functionB((x-(n+1)*h)/h)-4*functionB((x-(n+2)*h)/h)
	else:
		print("Valor de Matriz Invalido")

	return phi;

def functionPhiDerivada(x,j):
	if j==0:
		phi=functionBDerivada((x-0*h)/h)-4*functionBDerivada((x-(-1)*h)/h)
	elif j==1:
		phi=functionBDerivada((x-1*h)/h)-functionBDerivada((x-(-1)*h)/h)
	elif j>=2 and j<=n-1:
		phi=functionBDerivada((x-j*h)/h)
	elif j==n:
		phi=functionBDerivada((x-n*h)/h)-functionBDerivada((x-(n+2)*h)/h)
	elif j==n+1:
		phi=functionBDerivada((x-(n+1)*h)/h)-4*functionBDerivada((x-(n+2)*h)/h)
	else:
		phi=0;
		print("Valor de Matriz Invalido")
	return (phi);

n=9
h=1/(n+1)
vetorx = np.arange(h,1,h/10)
Plot=np.zeros(len(vetorx))
Plot1=np.zeros(len(vetorx))

for i in range (0,len(vetorx)):
	Plot[i]=functionPhi(vetorx[i],5)
	Plot1[i]=functionPhiDerivada(vetorx[i],5)


plt.plot(vetorx,Plot)
plt.plot(vetorx,Plot1)
plt.xlabel("Valores de X")
plt.ylabel("Função Phi")
plt.title("Função Phi")
plt.show()