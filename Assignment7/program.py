####################################################
####### Applied Programming - Assignment 7 #########
####### Done by Abishek S ##########################
####################################################

#Importing libraries
import pylab as pl
import numpy as np
import sys
import scipy.signal as sp
import matplotlib
from sympy import *

init_printing()

#Complex variable of Laplace equations
s = symbols('s')

#Defining a general utility function PLOT for plotting
def PLOT(x,y,fig_no = 0,label_x = r'$\rightarrow$',label_y = r'$\rightarrow$',fn = pl.plot,arg3 = 'b-',title = "Plot",grids = True,cmap = matplotlib.cm.jet,label = ''):
	'''Plotting function to make standard plots'''
	pl.figure(fig_no)
	pl.grid(grids)
	if fn == pl.contourf:
		fn(x,y,arg3,cmap = cmap)
		pl.colorbar()
	else:
		if label == '':
			fn(x,y,arg3)
		else:
			fn(x,y,arg3,label = label)
			pl.legend()
	pl.xlabel(label_x,size = 17)
	pl.ylabel(label_y,size = 17)
	pl.title(title)

#Defining a general utility function for making BODE Plots
def bodeplot(H,fig_no = 0,title = ''):
	'''Makes Bode Plots'''
	pl.figure(fig_no)
	pl.suptitle(title)
	w,s,phi = H.bode()
	pl.subplot(2,1,1)
	pl.semilogx(w,s)
	pl.xlabel(r'$\omega$',size=17)
	pl.ylabel(r'$|H(j\omega)|-(in dB)$',size =17)
	pl.subplot(2,1,2)
	pl.semilogx(w,phi)
	pl.xlabel(r'$\omega$',size=17)
	pl.ylabel(r'$\angle(H(j\omega))$',size =17)

#Solving the matrix equation and getting the output Voltage for Low Pass Filter
def LowPass(R1,R2,C1,C2,G,Vi=1):
    '''Active 2nd order low pass butterworth filter using opamp''' 
    A = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b = Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return(A,b,V)

#Solving the matrix equation and getting the output Voltage for High Pass Filter
def HighPass(R1,R2,C1,C2,G,Vi = 1):
    '''Active 2nd order high pass filter using opamp'''
    A = Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R2*C2)),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return(A,b,V)

#Extract the Transfer function and convert it to scipy lti form
def Simplify_H(V):
    '''Extracts Transfer function from the matrix inversion result'''
    Vo = V[3]  #The third element in the V column is the output voltage
    Vo = expand(simplify(Vo))  #converting to rational form
    H = SympyToScipy(Vo)
    return H

#Convert the Transfer function form sympy to scipy lti form
def SympyToScipy(Vo):
    '''Converts Transfer function in sympy to scipy'''
    v1 = fraction(Vo) #converting to numerator and denominator form
    n,d = Poly(v1[0],s),poly(v1[1],s)
    numer,denom = n.all_coeffs(), d.all_coeffs()  #extract the coefficients of 's'
    numer,denom = [float(f) for f in numer], [float(f) for f in denom]
    H = sp.lti(numer,denom)  #converting to scipy lti form
    return H

#Values of R1,R2,C1,C2,G chosen from given circuit
A,b,V = LowPass(10000,10000,1e-9,1e-9,1.586,1)
H_lp = Simplify_H(V)

#Values of R1,R2,C1,C2,G chosen from given circuit
A,b,V = HighPass(1e4,1e4,1e-9,1e-9,1.586,1)
H_hp = Simplify_H(V)

t = np.arange(0,1e-2,1e-7)  #Time scale

#Plotting Bode Plots of the high pass and low pass filters
bodeplot(H_lp,1,"Bode plot of Low Pass Filter")
pl.show()

bodeplot(H_hp,2,"Bode plot of High Pass Filter")
pl.show()


#Plotting step response of the two filters
A,b,V = LowPass(10000,10000,1e-9,1e-9,1.586,1/s)
_,vtd = sp.impulse(Simplify_H(V),None,t)
PLOT(t,vtd,3,r"$t$",r'$V_o(t)$',pl.plot,title = 'Step response of Low Pass Filter')
pl.show()

A,b,V = HighPass(1e4,1e4,1e-9,1e-9,1.586,1/s)
_,vtd = sp.impulse(Simplify_H(V),None,t)
PLOT(t,vtd,4,r"$t$",r'$V_o(t)$',pl.plot,title = 'Step response of High Pass Filter')
pl.show()

#Response of high and low pass filter to sum of sinusoids
#Frequency of one sinusoid is in the pass band of low pass filter and stop band of high pass filter
#Frequency of other sinusoid is in the pass band of high pass filter and stop band of low pass filter
t = np.arange(0,5e-3,1e-7)  #Time scale
inp = np.sin(2000*np.pi*t)+np.cos(2*(10**6)*np.pi*t)
PLOT(t,inp,5,r'$t$',r'$V_i(t)$',title = 'Sum of 2 sinusoids Input')
pl.show()

t,vtd,svec = sp.lsim(H_lp,inp,t)
PLOT(t,vtd,6,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of low pass filter to sum of 2 sinusoids')
pl.show()

t = np.arange(0,1e-5,1e-7)  #Time scale
inp = np.sin(2000*np.pi*t)+np.cos(2*(10**6)*np.pi*t)
t,vtd,svec = sp.lsim(H_hp,inp,t)
PLOT(t,vtd,7,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of high pass filter to sum of 2 sinusoids')
pl.show()

#Response of Low Pass and High Pass filter to high frequency damped sinusoid
damping_factor = 3000

t = np.arange(0,5e-4,1e-7)  #Time scale
inp = np.exp(-damping_factor*t)*(np.cos((10**6)*np.pi*t))
PLOT(t,inp,8,r'$t$',r'$V_i(t)$',title = 'High Frequency damped sinusoid input')
pl.show()

t,vtd,svec = sp.lsim(H_lp,inp,t)
PLOT(t,vtd,9,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of low pass filter to high frequency damped sinusoid')
pl.show()

t,vtd,svec = sp.lsim(H_hp,inp,t)
PLOT(t,vtd,10,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of high pass filter to high frequency damped sinusoid')
pl.show()

#Response of Low Pass and High Pass filter to low damped sinusoid
damping_factor = 100

t = np.arange(0,5e-2,1e-7)  #Time scale
inp = np.exp(-damping_factor*t)*(np.sin(1000*np.pi*t))
PLOT(t,inp,11,r'$t$',r'$V_i(t)$',title = 'Low Frequency damped sinusoid input')
pl.show()

t,vtd,svec = sp.lsim(H_lp,inp,t)
PLOT(t,vtd,12,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of low pass filter to low frequency damped sinusoid')
pl.show()

t,vtd,svec = sp.lsim(H_hp,inp,t)
PLOT(t,vtd,13,r"$t$",r'$V_o(t)$',pl.plot,title = 'Response of high pass filter to low frequency damped sinusoid')
pl.show()

############################ END OF PROGRAM ####################################