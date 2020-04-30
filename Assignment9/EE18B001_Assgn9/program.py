####################################################
####### Applied Programming - Assignment 9 #########
####### Done by Abishek S ##########################
####################################################

################################################################
#Importing Necessary Libraries 
from pylab import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

################################################################
#Defining utility functions for our ease

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
	#The sample corresponding to -tmax should be set zero
	y[0] = 0 
	Y = fftshift(fft(fftshift(y)))/float(N)
	w = linspace(-pi,pi,N,endpoint = False)
	#The range of frequencies
	w = w*(N/(end-st))

	fig, (ax1,ax2) = subplots(2,1)
	suptitle(f"Spectrum of {name}")
	ax1.plot(w,abs(Y),lw=1)
	ax1.set_ylabel(r"$|Y|$",size=16)
	ax1.grid(True)

	ax2.plot(w,angle(Y),'ro',lw=1)  
	ax2.grid(True)
	ax2.set_ylabel(r"Phase of $Y$",size=16)
	ax2.set_xlabel(r"$\omega$",size=16)

	return ax1,ax2,Y,w

################################################################
#Example 1 - sin(sqrt(2)*t)

y1 = lambda t: sin(sqrt(2)*t)
DFT(y1,(-pi,pi),64,r"$\sin(\sqrt{2}t)$")
#show()

####### Let's plot the extension and replication of our window

t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
# y=sin(sqrt(2)*t)
figure()
plot(t1,sin(sqrt(2)*t1),'b',lw=2)
plot(t2,sin(sqrt(2)*t2),'r',lw=2)
plot(t3,sin(sqrt(2)*t3),'r',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$")
grid(True)
#show()

y=sin(sqrt(2)*t1)
figure()
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)
#show()

####### Let's use windowing technique to get better peaks' amplitude

t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
y=t
#The sample corresponding to -tmax should be set zero
y[0]=0 
y=fftshift(y) #make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
semilogx(abs(w),20*log10(abs(Y)),lw=2)
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"],size=16)
ylabel(r"$|Y|$ (dB)",size=16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$",size=16)
grid(True)
#show()

wnd = lambda t: fftshift(0.54+0.46*cos(2*pi*t/len(t)))  #hamming window
t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
n=arange(64)
y=sin(sqrt(2)*t1)*wnd(n)
figure()
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$")
grid(True)
#show()

y2 = lambda t,n: sin(sqrt(2)*t)*wnd(arange(n))
y2_64 = lambda t : y2(t,64)
ax1,ax2,*_ = DFT(y2_64,(-pi,pi),64,r'$\sin(\sqrt{2}t*w(n)$)')
ax1.set_xlim(-8,8)
ax2.set_xlim(-8,8)
#show()

#Using 4 times the number of samples used earlier without affecting sampling frequency
y2_256 = lambda t : y2(t,256)
ax1,ax2,*_ = DFT(y2_256,(-4*pi,4*pi),256,r'$\sin(\sqrt{2}t*w(n))\ with\ better\ sampling$')
ax1.set_xlim(-4,4)
ax2.set_xlim(-4,4)
#show()

################################################################
#Question 2 - cos^3(wt)

#Without windowing
y3 = lambda t: (cos(0.86*t))**3
ax1,ax2,*_ = DFT(y3,(-pi,pi),64,r'$\cos^{3}(\omega_{0}t)$')
ax1.set_xlim(-10,10)
ax2.set_xlim(-10,10)
#show()

#With windowing
y3_w = lambda t: y3(t)*wnd(arange(256))
ax1,ax2,*_ = DFT(y3_w,(-4*pi,4*pi),256,r'$cos^{3}(\omega_{0}t)$')
ax1.set_xlim(-10,10)
ax2.set_xlim(-10,10)
#show()

################################################################
#Question 3 - arbitrary cos(wt + d)

def est_delta(w,Y,sup = 1e-3,window = 1):
	'''Estimates delta (d) from the spectrum of cos(w*t+d)'''
	ii_1 = where(logical_and(abs(Y)>sup, w>0))[0]
	sort(ii_1)
	points = ii_1[1:window+1]
	#weighted average for first 2 points
	return sum(angle(Y[points]))/len(points)

def est_omega(w,Y):
	'''Estimates omega (w) from the spectrum of cos(w*t+d)'''
	ii = where(w > 0)
	#omega estimated by weighted average
	return sum(abs(Y[ii])**2 * w[ii])/sum(abs(Y[ii])**2)

def CosEst(ww,d):
    fn = lambda t:  cos(ww*t+d)*wnd(arange(len(t)))
    ax1,ax2,Y,w = DFT(fn,(-pi,pi),128,r'$cos(\omega_{0}t+\delta)$')
    ax1.set_xlim(-10,10)
    ax2.set_xlim(-10,10)
    #show()

    print('Noiseless Signal parameters : ')
    print('\u03C9 :',est_omega(w,Y))
    print('\u03B4 :',est_delta(w,Y))
    return(Y)

Yf = CosEst(1.5,0.5)

################################################################
#Question 4 - the cos(wt + d) with gaussian noise added

def CosEst_Noisy(ww,d):
    fn = lambda t:  cos(ww*t+d)*wnd(arange(len(t))) + 0.1*randn(len(t))
    ax1,ax2,Y,w = DFT(fn,(-pi,pi),128,r'$cos(\omega_{0}t+\delta)\ with\ noise$')
    ax1.set_xlim(-10,10)
    ax2.set_xlim(-10,10)
    
    print('Noisy Signal parameters : ')
    print('\u03C9 :',est_omega(w,Y))
    print('\u03B4 :',est_delta(w,Y))
    return(Y)

Yf = CosEst_Noisy(1.5,0.5)
#show()

################################################################
#Question 5 - Chirped signal

#Plotting the time domain signal
ychirp = lambda t: cos(16*t*(1.5+t/(2*pi)))*wnd(arange(len(t)))
figure()
t = linspace(-pi,pi,1024)
plot(t,ychirp(t))
xlabel('t')
ylabel('y(t)')
title("Chirped Signal in time domain")
#show()

#Plotting the DFT spectrum
ax1,ax2,*_ = DFT(ychirp,(-pi,pi),1024,r'chirped signal')
ax1.set_xlim(-60,60)
ax2.set_xlim(-60,60)
#show()

################################################################
#Question 6 - Chirped signal - Time frequency plot

y_ch = cos(16*t*(1.5+t/(2*pi)))
NR = 64
NC = 1024//64
y_2D = zeros((NR,NC))
Y_2D = zeros((NR,NC),dtype = complex)
for i in range(NC):
	#Windowing
	y_2D[:,i] = y_ch[i*NR:(i+1)*NR]*wnd(arange(64))
	#The sample corresponding to -tmax should be set zero
	y_2D[0,i] = 0
	Y_2D[:,i] = fftshift(fft(fftshift(y_2D[:,i])))/float(NR)
x = linspace(-pi,pi,16,endpoint = False)
w = linspace(-pi,pi,64,endpoint = False)
w = w*1024.0/(2*pi);
#Forming the x and y values for the surface plot
wv,xv = meshgrid(w,x,indexing = 'ij')

#We plot the surface plot of magnitude of Y
fig = figure()
ax = p3.Axes3D(fig)
title('3-D spectrogram - Magnitude')
surf = ax.plot_surface(wv,xv,abs(Y_2D),cmap = cm.coolwarm,linewidth = 0,antialiased = False)
fig.colorbar(surf,shrink = 0.5,aspect = 5)
xlabel("Frequency")
ylabel("Time")
#show()

#We plot the contour plot of magnitude of Y
fig = figure()
title('Contour plot of Magnitude')
surf = contourf(xv,wv,abs(Y_2D))
ylim([-50,50])
ylabel("Frequency")
xlabel("Time")
fig.colorbar(surf)
show()

#We plot the surface plot of angle of Y
fig = figure()
ax = p3.Axes3D(fig)
title('3-D spectrogram - Angle')
surf = ax.plot_surface(wv,xv,angle(Y_2D),cmap = cm.coolwarm,linewidth = 0,antialiased = False)
fig.colorbar(surf,shrink = 0.5,aspect = 5)
xlabel("Frequency")
ylabel("Time")
#show()
############################ END OF PROGRAM ####################################