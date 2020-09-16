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
	return (phi/h);

def functionk(x):
   func = 1
   return func;

def functionq(x):
   func = 0
   return func;

def functionf(x):
   func = 12*x*(1-x)-2
   return func;

def function1(x,i,j):
  func= functionk(x)*functionPhiDerivada(x,i)*functionPhiDerivada(x,j)+functionq(x)*functionPhi(x,i)*functionPhi(x,j)
  return func;


def function2(x,i,j):
  func= functionf(x)*functionPhi(x,i)
  return func;

def trapez(function,a,b,t,i,posicaoi,posicaoj):
  #a e b dao o intervalo
  #t e o ultimo valor calculado
  #i e a quantidade de intervalos
  #newT e o valor do trapezio calculado
	h  = (b-a)/2**i
	if i==0:
	  #Se i for zero, calcula o primeiro trapezio
	  newT=(h/2)*(function(a,posicaoi,posicaoj)+function(b,posicaoi,posicaoj))
	else:     #Calculo dos trapezios, usando valor calculado previamente e pontos novos
		newPoints = np.arange(a+h,b,2*h)
		newFunctionsValues=np.zeros(len(newPoints))
		for i in range (0,len(newPoints)):
			newFunctionsValues[i]=function(newPoints[i],i=posicaoi,j=posicaoj)
		newTrapzSum=np.sum(newFunctionsValues)
		newT = (t/2) + h*newTrapzSum
	return newT


def romb(function,a,b,n,erro,ITMAX,posicaoi,posicaoj):

  matrixRomberg=np.zeros((n+1,n+1))
  iterations=0
  convergence=0

  for i in range(0,n+1):
    matrixRomberg[i][0]=trapez(function,a,b,matrixRomberg[i-1][0],i,posicaoi=posicaoi, posicaoj=posicaoj)
    for j in range (0,i):
      matrixRomberg[i][j+1] = 1.0/(4**(j+1)-1)*(4**(j+1)*matrixRomberg[i][j] - matrixRomberg[i-1][j])
    iterations=iterations+1
    if i>0:
      if abs(matrixRomberg[i][i]-matrixRomberg[i][i-1])<(erro*matrixRomberg[i][i]):
        convergence=1
        #print(matrixRomberg)
        return matrixRomberg[i][i],iterations,convergence
      elif iterations>=ITMAX:
        convergence=2
        #print(matrixRomberg)
        return matrixRomberg[i][i],iterations,convergence

  #print(matrixRomberg)
  return matrixRomberg[-1][-1],iterations,convergence

def resolucaodeMatrix(diagonalPrincipal,diagonal1,diagonal2,diagonal3,vetorb):
	s1Diagonal=np.diagflat(diagonal1,1)
	s2Diagonal=np.diagflat(diagonal2,2)
	s3Diagonal=np.diagflat(diagonal3,3)
	i1Diagonal=np.diagflat(diagonal1,-1)
	i2Diagonal=np.diagflat(diagonal2,-2)
	i3Diagonal=np.diagflat(diagonal3,-3)
	Diagonal=np.diagflat(diagonalPrincipal)

	matrixFinal=Diagonal+s1Diagonal+s2Diagonal+s3Diagonal+i1Diagonal+i2Diagonal+i3Diagonal

	matrixLD=np.zeros((n,n))
	matrixLD[0][0]=matrixFinal[0][0]

	for i in range (1,n):
		for j in range (0,i+1):
			if (i-j)>3:
				matrixLD[i][j]=0
			elif j==0:
				matrixLD[i][0]=matrixFinal[i][0]/matrixLD[0][0]
			elif i==j:
				somaTotal=0
				for k in range(0,j):
					somaTotal=somaTotal+(matrixLD[j][k]*matrixLD[j][k]*matrixLD[k][k])
				matrixLD[i][i]=matrixFinal[i][i]-somaTotal
			else:
				somaTotal=0
				for k in range(0,j):
					somaTotal=somaTotal+(matrixLD[i][k]*matrixLD[j][k]*matrixLD[k][k])
				matrixLD[i][j]=(1/matrixLD[j][j])*(matrixFinal[i][j]-somaTotal)

	matrixD=np.diagflat(matrixLD.diagonal())
	for i in range (0,n):
		matrixLD[i][i]=1
	LT=matrixLD.transpose()
	vetorz=np.zeros(n)
	vetory=np.zeros(n)
	vetorx=np.zeros(n)


	vetorz[0]=vetorb[0]

	for i in range (1,n):
		vetorz[i]=vetorb[i]-np.dot(matrixLD[i],vetorz)


	for i in range (0,n):
		vetory[i]=vetorz[i]/matrixD[i][i]


	vetorx[-1]=vetory[-1]

	for i in range (2,n+1):
		vetorx[-i]=vetory[-i]-np.dot(LT[-i],vetorx)

	return vetorx

def functionU(alphai,posicaoi,x):
	func = alphai*functionPhi(x,posicaoi)
	return func

def calculoAlpha(n):


	vetorb=np.zeros(n)
	diagonalPrincipal=np.zeros(n)
	diagonal1=np.zeros(n-1)
	diagonal2=np.zeros(n-2)
	diagonal3=np.zeros(n-3)


	for i in range (0,n):
		U=min(0,(i+3)*h)
		L=max(1,(i-1)*h)
		diagonalPrincipal[i],iteracoes,convergencia=romb(function=function1,a=L,b=U,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1),posicaoj=(i+1))
		vetorb[i],iteracoes,convergencia=romb(function=function2,a=L,b=U,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+1),posicaoj=(i+1))


	for i in range (0,n-1):
		U=min(0,(i+4)*h)
		L=max(1,(i-1)*h)
		diagonal1[i],iteracoes,convergencia=romb(function=function1,a=L,b=U,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+2),posicaoj=(i+1))

	for i in range (0,n-2):
		U=min(0,(i+5)*h)
		L=max(1,(i-1)*h)
		diagonal2[i],iteracoes,convergencia=romb(function=function1,a=L,b=U,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+3),posicaoj=(i+1))

	for i in range (0,n-3):
		U=min(0,(i+5)*h)
		L=max(1,(i-1)*h)
		diagonal3[i],iteracoes,convergencia=romb(function=function1,a=L,b=U,n=n,erro=1e-10,ITMAX=10,posicaoi=(i+4),posicaoj=(i+1))

	vetoralpha=resolucaodeMatrix(diagonalPrincipal,diagonal1,diagonal2,diagonal3,vetorb)
	#print (vetoralpha) #- caso queira ver o vetor alfa
	return vetoralpha

def calculoErro(n):
	valoresdealpha=calculoAlpha(n)
	vetorU=np.zeros(10*n)
	vetorProvisorio=np.zeros(n)
	vetorFuncao=np.zeros(10*n)
	vetorErro=np.zeros(10*n)

	for j in range (0,10*n):
		for i in range (0,n):
			vetorProvisorio[i]=functionU(valoresdealpha[i],(i+1),(j+1)*h/10)
		vetorU[j]=np.sum(vetorProvisorio)



	for i in range (0,10*n):
		vetorFuncao[i]=((i+1)/10*h)**2*(((i+1)/10*h)-1)**2
		vetorErro[i]=abs(vetorFuncao[i]-vetorU[i])

	return vetorErro.max()

valoresDeN=np.array([15,31,63,127,255])
valoresDeErro=np.zeros(len(valoresDeN))
for i in range (0,len(valoresDeN)):
  h=1/(valoresDeN[i]+1)
  n=valoresDeN[i]
  valoresDeErro[i]=calculoErro(valoresDeN[i])
  print(" Para o valor de n: ",valoresDeN[i], "temos o erro: ", valoresDeErro[i])

plt.yscale('log',basey=10) 
plt.plot(valoresDeN,valoresDeErro)
plt.xlabel("Valores de N")
plt.ylabel("Erro")
plt.title("Erro por N")
plt.show()