####################################################
####### Applied Programming - Assignment 6 #########
####### Done by Abishek S ##########################
####################################################

#Importing libraries
import pylab as pl
import numpy as np
import sys
import scipy.signal as sp
import matplotlib

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
def bodeplot(w,s,phi):
	'''Makes Bode Plots'''
	pl.subplot(2,1,1)
	pl.semilogx(w,s)
	pl.xlabel(r'$\omega$',size=17)
	pl.ylabel(r'$|H(j\omega)|-(in dB)$',size =17)
	pl.subplot(2,1,2)
	pl.semilogx(w,phi)
	pl.xlabel(r'$\omega$',size=17)
	pl.ylabel(r'$\angle(H(j\omega))$',size =17)

################################################################
#Qns1,2 - Solving for x for given f(t) and differential equation

def input_fn(decay=0.5,cos_term=1.5):
	'''Laplace Tranform of input function'''
	return (np.poly1d([1,decay]),np.poly1d([1,2*decay,decay**2 + cos_term**2]))

def transfer_fn(a=1,b=0,c=2.25):
	'''Transfer function of system
	for (a*s^2 + b*s + c)*X(s)'''
	return (np.poly1d([1]),np.poly1d([a,b,c]))

def zero_st(a=1,b=0,c=2.25,xi=0,xii=0):
	'''Zero state Response of system
	for (a*s^2 + b*s + c)*X(s) with x(0) = xi and dx/dt at t = 0 is xii'''
	return (np.poly1d([a*xi,(b*xi)+(a*xii)]),np.poly1d([a,b,c]))

def output_fn(F,a=1,b=0,c=2.25,xi=0,xii=0):
	'''Laplace Transform of system output
	where F is input in Laplace domain'''
	ziN,ziD = transfer_fn(a,b,c)
	zsN,zsD = zero_st(a,b,c,xi,xii)
	return (np.polyadd(np.polymul(F[0],ziN),np.polymul(F[1],zsN)),np.polymul(F[1],ziD))

t = np.linspace(0,50,200)

#Plotting the output (x) for decay rate = 0.5 (Qn1)
t,x1 = sp.impulse(output_fn(input_fn()),None,t)
PLOT(t,x1,1,r'$t\rightarrow$',r'$x(t)\rightarrow$',title = 'System response for decay rate = 0.5')
pl.show()

#Plotting the output (x) for decay rate = 0.05 (Qn2)
t,x2 = sp.impulse(output_fn(input_fn(decay=0.05)),None,t)
PLOT(t,x2,2,r'$t\rightarrow$',r'$x(t)\rightarrow$',title = 'System response for decay rate = 0.05')
pl.show()

################################################################
#Qn3 - Deriving Transfer function and using convolutiom to find the system output

def input_td(t,decay=0.5,cos_term=1.5):
	'''Return list of input values in time domain corresponding to the t values given'''
	return (np.cos(f*t)*np.exp(-1*decay*t))

t = np.linspace(0,100,300)
H = transfer_fn()

#We use the following 5 colours for the lines :
#Black Green Red Cyan and Magenta
for f,col in zip(np.arange(1.4,1.61,0.05),['k','g','r','c','m']):
	u = input_td(t,0.05,f)
	t,y,svec = sp.lsim(H,u,t)
	PLOT(t,y,3,r'$t$',r'$x(t)$',title = 'System responses for different frequencies with decay rate = 0.05',arg3 = col+'-',label = 'freq '+str(f))

pl.show()

#We make the Bode plot of the transfer function to understand it better
w,s,phi = sp.lti(H[0],H[1]).bode()
pl.figure(4)
bodeplot(w,s,phi) 
pl.suptitle('Bode plot of Transfer function of single spring system')
pl.show()

################################################################
#Qn4 - Solving coupled spring system

t = np.linspace(0,20,200)
#Transfer function for x
X = sp.lti([1,0,2],[1,0,3,0])

#Transfer function for y
Y = sp.lti([2],[1,0,3,0])

t,x = sp.impulse(X,None,t)
t,y = sp.impulse(Y,None,t)

#Plotting x(t) and y(t) in a single graph
PLOT(t,x,5,'t','f(t)',title = 'Plot of x(t) and y(t) in coupled spring system',label = 'x(t)')
PLOT(t,y,5,'t','f(t)',title = 'Plot of x(t) and y(t) in coupled spring system',arg3 = 'g-',label = 'y(t)')
pl.show()

################################################################
#Qns5,6 - Solving RLC circuit which is analogous to the spring system

def RLCtf():
	'''Transfer function of the RLC circuit given'''
	return (np.poly1d([1]),np.poly1d([1e-12,1e-4,1]))

def RLCinp():
	'''Input to the RLC circuit'''
	return (np.cos(1e3*t)-np.cos(1e6*t))

#Defining time vector appropriately from 0 to 10msec with 1e-7 time steps to cpature the fast variation
t = np.arange(0,10e-3,1e-7)
H = RLCtf()
#We make the Bode plot of the transfer function
w,s,phi = sp.lti(H[0],H[1]).bode()
pl.figure(4)
bodeplot(w,s,phi) 
pl.suptitle('Bode plot of Transfer function of RLC filter')
pl.show()

#Solving for the RLC filter output
u = RLCinp()
t,x,svec = sp.lsim(H,u,t)

#Plotting the slow time and fast time outputs separately
PLOT(t,x,6,'t',r'$V_{o}(t)$',title = 'Plot of Output voltage of RLC circuit - Slow time')
pl.show()
PLOT(t[:300],x[:300],7,'t',r'$V_{o}(t)$',title = 'Plot of Output voltage of RLC circuit - Fast time')
pl.show()

########################### END OF PROGRAM #####################################