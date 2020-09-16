import numpy as np
from matplotlib import pyplot as plt


def functionPhi(x,i):
	if x>=(i-1)*h and x<=(i*h):
		phi=(x-(i-1)*h)/h
	elif x>i*h and x<=(i+1)*h:
		phi=((i+1)*h-x)/h
	else:
		phi=0
	return phi


def functionPhiDerivada(x,i):
  if x>=(i-1)*h and x<=(i*h):
  	phi=1/h
  elif x>i*h and x<=(i+1)*h:
  	phi=-1/h
  else:
  	phi=0
  return phi

def functionk(x):
   sigma=0.1
   func = 2*((np.exp(1))**(((-1*(x-ComprimentoL/2)**2))/(sigma)**2))
   return func;

def functionkderivada(x):
   sigma=0.1
   func = 2*((np.exp(1))**(((-1*(x-ComprimentoL/2)**2))/(sigma)**2))*((2*x-ComprimentoL)/(sigma)**2)
   return func;

def functionq(x):
   func = 0
   return func;

def functionf(x):
  f0=1
  d=0.4
  if x>((ComprimentoL/2)-d) and x<((ComprimentoL/2)+d):
    fx=f0
  else:
    fx=0
  func= fx + functionkderivada(x)*(valorfinal-valorinicial)-functionq(x)*(valorinicial+(valorfinal-valorinicial)*x)
  return func;

def functionQ1(x,i):
  func= functionq(x)*((i+1)*h-x)*(x-i*h)
  return func;

def functionQ2(x,i):
  func= functionq(x)*(x-(i-1)*h)**2
  return func;

def functionQ3(x,i):
  func= functionq(x)*(((i+1)*h)-x)**2
  return func;

def functionQ4(x,i):
  func= functionk(x)
  return func;

def functionQ5(x,i):
  func= functionf(x)*(x-(i-1)*h)
  return func;

def functionQ6(x,i):
  func= functionf(x)*(((i+1)*h)-x)
  return func;

def trapez(function,a,b,t,i,posicaoi):
  #a e b dao o intervalo
  #t e o ultimo valor calculado
  #i e a quantidade de intervalos
  #newT e o valor do trapezio calculado
  h  = (b-a)/2**i
  if i==0:
    #Se i for zero, calcula o primeiro trapezio
    newT=(h/2)*(function(a,posicaoi)+function(b,posicaoi))
  else:     #Calculo dos trapezios, usando valor calculado previamente e pontos novos
    newPoints = np.arange(a+h,b,2*h)
    newFunctionsValues=np.zeros(len(newPoints))
    for i in range (0,len(newPoints)):
  	   newFunctionsValues[i]=function(newPoints[i],i=posicaoi)
    newTrapzSum=np.sum(newFunctionsValues)
    newT = (t/2) + h*newTrapzSum
  return newT


def romb(function,a,b,n,erro,ITMAX,posicaoi):

  matrixRomberg=np.zeros((n+1,n+1))
  iterations=0
  convergence=0

  for i in range(0,n+1):
    matrixRomberg[i][0]=trapez(function,a,b,matrixRomberg[i-1][0],i,posicaoi=posicaoi)
    for j in range (0,i):
      matrixRomberg[i][j+1] = 1.0/(4**(j+1)-1)*(4**(j+1)*matrixRomberg[i][j] - matrixRomberg[i-1][j])
    iterations=iterations+1
    if i>0:
      if abs(matrixRomberg[i][i]-matrixRomberg[i][i-1])<(erro*matrixRomberg[i][i]):
        convergence=1
        #print(matrixRomberg[:,0])
        return matrixRomberg[i][i],iterations,convergence
      elif iterations>=ITMAX:
        convergence=2
        #print(matrixRomberg[:,0])
        return matrixRomberg[i][i],iterations,convergence

  #print(matrixRomberg[:,0])
  return matrixRomberg[-1][-1],iterations,convergence

def Thomas (a,b,d):
  z=len(a)
  for i in range (1,z):
    #Calcula a diagonal principal e o novo d
      a[i]=a[i]-(b[i-1]*b[i-1]/a[i-1])
      d[i]=d[i]-(b[i-1]*d[i-1]/a[i-1])
  
  #Calculo do vetor x
  x=np.zeros(z)
  x[-1]=d[-1]/a[-1]

  for i in range (z-2,-1,-1):
    x[i]=(d[i]-(b[i]*x[i+1]))/a[i]
  
  return x

def calculoDeX(n):

  tolerancia=1e-10

  Q1=np.zeros(n)
  Q2=np.zeros(n)
  Q3=np.zeros(n)
  Q4=np.zeros(n+1)
  Q5=np.zeros(n)
  Q6=np.zeros(n)
  diagonalPrincipal=np.zeros(n)
  diagonal1=np.zeros(n-1)
  vetorb=np.zeros(n)

  for i in range(0,n-1):
    rombergQ1,iteracoes,convergenciaQ1=romb(function=functionQ1,a=(i+1)*h,b=(i+2)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q1[i]=((1/h)**2)*rombergQ1

  for i in range(0,n):
    rombergQ2,iteracoes,convergenciaQ2=romb(function=functionQ2,a=(i)*h,b=(i+1)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q2[i]=((1/h)**2)*rombergQ2
    rombergQ3,iteracoes,convergenciaQ3=romb(function=functionQ3,a=(i+1)*h,b=(i+2)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q3[i]=((1/h)**2)*rombergQ3
    rombergQ5,iteracoes,convergenciaQ5=romb(function=functionQ5,a=(i)*h,b=(i+1)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q5[i]=(1/h)*rombergQ5
    rombergQ6,iteracoes,convergenciaQ6=romb(function=functionQ6,a=(i+1)*h,b=(i+2)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q6[i]=(1/h)*rombergQ6

  for i in range(0,n+1):
    rombergQ4,iteracoes,convergenciaQ4=romb(function=functionQ4,a=(i)*h,b=(i+1)*h,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1))
    Q4[i]=((1/h)**2)*rombergQ4

  for i in range(0,n-1):
    diagonal1[i]=Q1[i]-Q4[i+1]
    diagonalPrincipal[i]=Q4[i]+Q4[i+1]+Q2[i]+Q3[i]
    vetorb[i]=Q5[i]+Q6[i]


  diagonalPrincipal[n-1]=Q4[n-1]+Q4[n]+Q2[n-1]+Q3[n-1]
  vetorb[n-1]=Q5[n-1]+Q6[n-1]



  resultado=Thomas(diagonalPrincipal,diagonal1,vetorb)

  return resultado

valorinicial=0
valorfinal=0 #kappa

n=99
h=1/(n+1)
ComprimentoL=1

vetorv=calculoDeX(n)
vetoru=np.zeros(n)
valoresdex=np.arange(h,1,h)

for i in range (0,n):
  vetoru[i]=vetorv[i]+valorinicial+(valorfinal-valorinicial)*(i+1)*h
 

#print(vetoru)

plt.plot(valoresdex,vetoru)
plt.xlabel("Valores de X")
plt.ylabel("Dispersao de Poluente")
plt.title("Dispersao de Poluente")
plt.show()