####################################################
####### Applied Programming - Assignment 8 #########
####### Done by Abishek S ##########################
####################################################

################################################################
#Importing Necessary Libraries 
from pylab import *
import sys
import matplotlib

################################################################
#Example 1 - FFT and IFFT Basic

x=rand(100)
X=fft(x)
y=ifft(X)
print(c_[x,y][:10])  #printing first 10 lines alone to get an idea of how close x and y values are
print('Maximum absolute error between x and y is :',abs(x-y).max())

################################################################
#Example 2 - FFT of sin(5x) - Trial 1

x=linspace(0,2*pi,128)
y=sin(10*x)
Y=fft(y)
figure()
subplot(2,1,1)
plot(abs(Y),lw=2)  #line width of 2
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin(5t)$")
grid(True)  #for the grids to be shown
subplot(2,1,2)
plot(unwrap(angle(Y)),lw=2)
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
#savefig("ex2_i.png")  #Automatically saving the figure
show()

################################################################
#Example 2 - FFT of sin(5x) - Trial 2
#Since we didn't get the exact digital spectrum we expected

x=linspace(0,2*pi,129);x=x[:-1]  #last point is excluded
y=sin(5*x)
Y=fftshift(fft(y))/128.0  #fftshift converts from [0,2pi] to [-pi,pi] 
w=linspace(-64,63,128)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin(5t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)  #highlighting points where the magnitude is significant
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
#savefig("ex2_f.png")
show()

################################################################
#Example 3 - AM modulation of given function - Trial 1

t=linspace(0,2*pi,129);t=t[:-1]  #Low number of samples
y=(1+0.1*cos(t))*cos(10*t)  #AM with carrier at 10 and modulating freq of 1
Y=fftshift(fft(y))/128.0
w=linspace(-64,63,128)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("eg3_i.png")
show()

################################################################
#Example 3 - AM modulation of given function - Trial 2

t=linspace(-4*pi,4*pi,513);t=t[:-1]  #High number of samples, hence tighter spectrum
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-15,15])
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("eg3_f.png")
show()

################################################################
#Questions 2,3 - Obtaining spectrum for specified signals

def DFT(y_fn,tim = (-4*pi,4*pi),N = 512,name = ''):
	''' Utility function for generating spectrum of given function.

	y_fn : function for which spectrum needs to be obtained
	tim : time interval of form (starting time,ending time) where ending time has to be excluded
	N : number of samples in time domain
	name : Name of function for the spectrum
	'''
	st,end = tim

	t = linspace(st,end,N,endpoint = False)
	y = y_fn(t)
	Y = fftshift(fft(y))/float(N)
	w = linspace(-pi,pi,N,endpoint = False)
	w = w*(N/(end-st))
	fig, (ax1,ax2) = subplots(2,1)
	ii = where(abs(Y)>1e-3)
	ax1.plot(w,abs(Y),lw=1)
	ax1.set_xlim([-2*max(w[ii]),2*max(w[ii])])
	ax1.set_ylabel(r"$|Y|$",size=16)
	suptitle(f"Spectrum of {name}")
	ax1.grid(True)
	#ax2.plot(w,angle(Y),'ro',lw=1)
	ax2.plot(w[ii],angle(Y[ii]),'go',lw=1)  #Plotting only the phase of relevant points
	ax2.grid(True)
	ax2.set_xlim([-2*max(w[ii]),2*max(w[ii])])
	ax2.set_ylabel(r"Phase of $Y$",size=16)
	ax2.set_xlabel(r"$\omega$",size=16)
	show()

y1 = lambda t : (sin(t))**3
y2 = lambda t : (cos(t))**3
y3 = lambda t : cos(20*t+5*cos(t))

DFT(y1,(-4*pi,4*pi),256,r'$sin^{3}t$')
DFT(y2,(-2*pi,2*pi),256,r"$cos^{3}t$")
DFT(y3,(-4*pi,4*pi),2048,r"$cos(20t+5cos(t))$")

################################################################
#Question 4 - Obtaining spectrum for Gaussian (which is not Bandlimited)

def estCTFT(y_fn,org_fn,tim = (-4*pi,4*pi),N = 512,name = ''):
    ''' Utility function to plot DFT for Gaussian by taking a window of values and estimate how close the DFT is to the CTFT.

	y_fn : function for which spectrum needs to be obtained
	org_fn : the orginal CTFT function
	tim : time interval of form (starting time,ending time) where ending time has to be excluded
	N : number of samples in time domain
	name : Name of function for the spectrum
	'''
    start,end = tim

    t = linspace(end,start,N,endpoint=False)
    y = y_fn(t)
    Y = fftshift(fft(ifftshift(y)))*(end-start)/(2*pi*N)  #Using ifftshift to remove phase issues

    w=linspace(-pi,pi,N,endpoint= False);
    w = w*(N/(end-start))
    
    ctft = org_fn(w)
    #Sum of absolute error between orginal CTFT and the DFT
    error = sum(abs(ctft-Y))
    print('Total absolute error between CTFT and DFT is :',error)
    
    fig, (ax1, ax2) = subplots(2, 1)
    ii = where(abs(Y) > 1e-3)
    suptitle("Spectrum of {}".format(name))
    #Plotting the orginal CTFT of Gaussian in the same plot
    ax1.plot(w,abs(ctft),'y-',label = 'CTFT')
    ax2.plot(w,angle(ctft),'yo',label = 'CTFT')

    ax1.plot(w,abs(Y),lw=1,label = 'DFT')
    ax1.set_xlim([-2*max(w[ii]),2*max(w[ii])])
    ax1.set_ylabel(r"$|Y|$",size=16)
    ax1.grid(True)
    ax2.plot(w[ii],angle(Y[ii]),'ro',label = 'DFT')  #Plotting only the phase of relevant points
    ax2.set_xlim([-2*max(w[ii]),2*max(w[ii])])  
    ax2.set_ylabel(r"Phase of $Y$",size=16)
    ax2.set_xlabel(r"$\omega$",size=16)
    ax2.grid(True)
    ax1.legend(loc = 'upper right')
    ax2.legend(loc = 'upper right')

    show()

y_gaussian = lambda t : exp(t**2/-2)
ctft_gaussian = lambda w : exp(-w**2/2)/sqrt(2*pi)

estCTFT(y_gaussian,ctft_gaussian,(-4*pi,4*pi),1024,r"$exp(-t^{2}/2)$")


############################ END OF PROGRAM ####################################