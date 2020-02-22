####################################################
####### Applied Programming - Assignment 3 #########
####### Done by Abishek S ##########################
####################################################

#Importing libraries
from scipy.integrate import quad
from scipy.linalg import lstsq
import pylab as pl
import numpy as np
import sys
import math

################################################################
#Part1 - Defining two functions and plotting them over [-2Pi,4Pi)

#Defining original Functions
def exp(x):
	return np.exp(x) if np.isscalar(x) else np.exp(np.array(x))

def coscos(x):
	return np.cos(np.cos(x)) if np.isscalar(x) else np.cos(np.cos(np.array(x)))

#Defining two functions' Periodic extension
def f1(x):
	X = np.asarray([x]) if np.isscalar(x) else np.asarray(x)
	X = list(map(lambda y: y - math.floor(y/(2*np.pi))*(2*np.pi),X))
	return np.exp(X[0]) if np.isscalar(x) else np.exp(np.array(X))

def f2(x):
	X = np.asarray([x]) if np.isscalar(x) else np.asarray(x)
	X = list(map(lambda y: y - math.floor(y/2*np.pi)*(2*np.pi),X))
	return np.cos(np.cos(X[0])) if np.isscalar(x) else np.cos(np.cos(np.array(X)))

X_list = pl.linspace(-2*np.pi,4*np.pi,300,endpoint = True)

#Plotting exp(x) original and periodic extensions in same graph
pl.figure(1)
pl.title('Original and periodic extension of exp(x)')
pl.semilogy(X_list,exp(X_list),label = 'exp(x) - Original')
pl.semilogy(X_list,f1(X_list),label = 'exp(x) - Periodic Extension')
pl.xlabel(r'X $\rightarrow$')
pl.ylabel(r'Y $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

#Plotting exp(x) original and periodic extensions in same graph
pl.figure(2)
pl.title('Original and periodic extension of cos(cos(x))')
pl.plot(X_list,coscos(X_list),label = 'cos(cos(x)) - Original')
pl.plot(X_list,f2(X_list),label = 'cos(cos(x)) - Periodic Extension')
pl.xlabel(r'X $\rightarrow$')
pl.ylabel(r'Y $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

################################################################
#Part2 - Integration using quad

def u(x,k,f):
	return f(x)*np.cos(k*x)

def v(x,k,f):
	return f(x)*np.sin(k*x)

a = np.zeros((2,26))
b = np.zeros((2,26))

for i in range(26):
	a[0][i] = quad(u,0,2*np.pi,args=(i,f1))[0]/np.pi
	if(i == 0):
		a[0][i] /= 2

for i in range(26):
	b[0][i] = quad(v,0,2*np.pi,args=(i,f1))[0]/np.pi
	if(i == 0):
		b[0][i] /= 2

#print('Printing Fourier series coefficients of exp(x)')
#print(a[0],end = '\n\n\n')
#print(b[0],end = '\n\n\n')

for i in range(26):
	a[1][i] = quad(u,0,2*np.pi,args=(i,f2))[0]/np.pi
	if(i == 0):
		a[1][i] /= 2

for i in range(26):
	b[1][i] = quad(v,0,2*np.pi,args=(i,f2))[0]/np.pi
	if(i == 0):
		b[1][i] /= 2


#print('Printing Fourier series coefficients of cos(cos(x))')
#print(a[1],end = '\n\n\n')
#print(b[1],end = '\n\n\n')

################################################################
#Part3 - Plotting semilogy and loglog plots of coefficients

#Forming the coefficients matrix
C = np.zeros((2,51))
C[0][0],C[1][0] = a[0][0],a[1][0]
for i in range(1,26):
	C[0][2*i-1],C[0][2*i] = a[0][i],b[0][i]
	C[1][2*i-1],C[1][2*i] = a[1][i],b[1][i]

#Plotting
pl.figure(3)
pl.title('Loglog Fourier coefficients of exp(x)')
pl.loglog(np.arange(1,52,1),np.absolute(C[0]),'ro',label = 'exp(x)')
pl.xlabel(r'log n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(4)
pl.title('Semilogy Fourier coefficients of exp(x)')
pl.semilogy(np.arange(1,52,1),np.absolute(C[0]),'ro',label = 'exp(x)')
pl.xlabel(r'n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(5)
pl.title('Loglog Fourier coefficients of cos(cos(x))')
pl.loglog(np.arange(1,52,1),np.absolute(C[1]),'ro',label = 'cos(cos(x))')
pl.xlabel(r'log n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(6)
pl.title('Semilogy Fourier coefficients of cos(cos(x))')
pl.semilogy(np.arange(1,52,1),np.absolute(C[1]),'ro',label = 'cos(cos(x))')
pl.xlabel(r'n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()


################################################################
#Parts 4,5 - Using lstsq to predict the Fourier series coefficients

x=pl.linspace(0,2*np.pi,400,endpoint = True)
# f1,f2 have been written to take a vector 
b1=f1(x) 
b2=f2(x)
#Providing the endpoint for f1(x) alone to make sure fourier series
#converges to the midpoint of left and right limits at a discontinuity
b1[-1] = exp(2*np.pi)
# allocate space for A
A=np.zeros((400,51))
#col1 is all ones
A[:,0]=1
for k in range(1,26):
	A[:,2*k-1]=np.cos(k*x) # cos(kx) column
	A[:,2*k]=np.sin(k*x)   # sin(kx) colum
#endfor

c1=lstsq(A,b1)[0]    # the ’[0]’ is to pull out the best fit vector
#lstsq returns a list
c2=lstsq(A,b2)[0]

#Plotting
pl.figure(3)
pl.title('Loglog Fourier coefficients of exp(x)')
pl.loglog(np.arange(1,52,1),np.absolute(C[0]),'ro',label = 'True')
pl.loglog(np.arange(1,52,1),np.absolute(c1),'go',label = 'Predicted')
pl.xlabel(r'log n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(4)
pl.title('Semilogy Fourier coefficients of exp(x)')
pl.semilogy(np.arange(1,52,1),np.absolute(C[0]),'ro',label = 'True')
pl.semilogy(np.arange(1,52,1),np.absolute(c1),'go',label = 'Predicted')
pl.xlabel(r'n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(5)
pl.title('Loglog Fourier coefficients of cos(cos(x))')
pl.loglog(np.arange(1,52,1),np.absolute(C[1]),'ro',label = 'True')
pl.loglog(np.arange(1,52,1),np.absolute(c2),'go',label = 'Predicted')
pl.xlabel(r'log n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

pl.figure(6)
pl.title('Semilogy Fourier coefficients of cos(cos(x))')
pl.semilogy(np.arange(1,52,1),np.absolute(C[1]),'ro',label = 'True')
pl.semilogy(np.arange(1,52,1),np.absolute(c2),'go',label = 'Predicted')
pl.xlabel(r'n $\rightarrow$')
pl.ylabel(r'log of Magnitude of coefficient $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

################################################################
#Part6 - Finding largest deviation between true and predicted coefficients

e1 = np.absolute(C[0]-c1)
e2 = np.absolute(C[1]-c2)

print("The maximum deviation for exp(x) function :",np.amax(e1))
print('The coefficients with largest deviation is/are :')

for i in np.argwhere(e1 == np.amax(e1)) : #Returning all maximum deviations
	coef = ''
	if i == 0:
		print('A0')
		continue
	coef += chr(65 + 1-(i%2))
	coef += str(int((i+1)/2))
	print(coef)

print("The maximum deviation for cos(cos(x)) function :",np.amax(e2))
print('The coefficients with largest deviation is/are :')

for i in np.argwhere(e1 == np.amax(e1)) : #Returning all maximum deviations
	coef = ''
	if i == 0:
		print('A0')
		continue
	coef += chr(65 + 1-(i%2))
	coef += str(int((i+1)/2))
	print(coef)

################################################################
#Part7 - Finding largest deviation between true and predicted coefficients

Acexp = np.matmul(A,c1)
Accos = np.matmul(A,c2)

new_x = pl.linspace(0,2*np.pi,400,endpoint = True)

#Plotting exp(x) predicted from lstsq in Fig 1
pl.figure(1)
pl.title('Original and lstsq predicted versions of exp(x)')
pl.semilogy(new_x,Acexp,'go',label = 'exp(x) - Predicted')
pl.semilogy(X_list,exp(X_list),label = 'exp(x) - Original')
pl.semilogy(X_list,f1(X_list),label = 'exp(x) - Periodic Extension')
pl.xlabel(r'X $\rightarrow$')
pl.ylabel(r'Y $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

#Plotting cos(cos(x)) predicted from lstsq in Fig 2
pl.figure(2)
pl.title('Original and lstsq predicted versions of cos(cos(x))')
pl.plot(new_x,Accos,'go',label = 'cos(cos(x)) - Predicted')
pl.plot(X_list,coscos(X_list),label = 'cos(cos(x)) - Original')
pl.plot(X_list,f2(X_list),label = 'cos(cos(x)) - Periodic Extension')
pl.xlabel(r'X $\rightarrow$')
pl.ylabel(r'Y $\rightarrow$')
pl.legend()
pl.grid(True)
pl.show()

############################# END OF PROGRAM ###################################
