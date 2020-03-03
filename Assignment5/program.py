####################################################
####### Applied Programming - Assignment 3 #########
####### Done by Abishek S ##########################
####################################################

#Importing libraries
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
import argparse
import numpy as np
from scipy.linalg import lstsq
import sys
import matplotlib

################################################################
#Setting Parameters

#Using argparse to create Nx,Ny,radius,Ni as optional command line arguments
parser = argparse.ArgumentParser(description='Defining Parameters')
parser.add_argument('-Nx',dest='Nx',type=int,default=25,help = 'Size along x')
parser.add_argument('-Ny',dest='Ny',type=int,default=25,help = 'Size along y')
parser.add_argument('-r',dest='radius',type=int,default=8,help = 'Radius of central lead')
parser.add_argument('-Ni',dest='Niter',type=int,default=1500,help = 'Number of interations to perform')
args = parser.parse_args()
Nx,Ny,Niter,radius = args.Nx,args.Ny,args.Niter,args.radius

################################################################
#Allocate Potential array and initialize it

#Phi is the potential array
Phi = np.zeros(shape = (Ny,Nx))
#Getting the x and y co-ordinates of the potential array
x = np.linspace(-0.5,0.5,Nx)
y = np.linspace(-0.5,0.5,Ny)
Y,X = np.meshgrid(y,x)

#Defining a general utility function PLOT for plotting
def PLOT(x,y,fig_no = 0,label_x = r'$\rightarrow$',label_y = r'$\rightarrow$',fn = pl.plot,arg3 = 'b-',title = "Plot",grids = True,cmap = matplotlib.cm.jet,label = ''):
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
	pl.xlabel(label_x)
	pl.ylabel(label_y)
	pl.title(title)

ii = np.where(X*X + Y*Y <= 0.35*0.35)
Phi[ii] = 1

#Plotting the contour of potential
PLOT(x,y,1,"X-axis","Y-axis",pl.contourf,Phi,"Contour Plot of Potential",True,matplotlib.cm.hot)
pl.plot(X[X*X + Y*Y < 0.35*0.35],Y[X*X + Y*Y < 0.35*0.35],'ro',label = 'Points with potential 1V')
pl.legend()
pl.show()

################################################################
#Iterating many times for convergence of Potential

errors = np.zeros(Niter)
for k in range(Niter):
	oldphi = Phi.copy()
	#Interior Points - average of surrounding points
	Phi[1:-1,1:-1] = 0.25*(Phi[1:-1,0:-2]+Phi[1:-1,2:]+Phi[0:-2,1:-1]+Phi[2:,1:-1])
	Phi[0,:] = 0  #Boundary condition for bottom side
	Phi[1:-1,0] = Phi[1:-1,1]  #Boundary condition for left side
	Phi[1:-1,-1] = Phi[1:-1,-2]  #Boundary condition for right side
	Phi[-1,:] = Phi[-2,:]  #Boundary condition for top side
	Phi[ii] = 1.0  #Condition for 1V points
	errors[k] = (abs(Phi-oldphi)).max()

#Plotting Errors
#Semilog Plot
PLOT(np.arange(1,Niter+1),errors,2,"Iteration Number","Log of Error",fn = pl.semilogy,title='Semilogy Plot of Error vs Iteration Number',label = 'Line')
pl.semilogy(np.arange(1,Niter+1,30),errors[0:Niter:30],'yo',label = 'Points')
pl.legend()
pl.show()

#Loglog plot
PLOT(np.arange(1,Niter+1),errors,3,"Log of Iteration Number","Log of Error",fn = pl.loglog,title='Loglog Plot of Error vs Iteration Number',label = 'Line')
pl.loglog(np.arange(1,Niter+1,30),errors[0:Niter:30],'yo',label = 'Points')
pl.legend()
pl.show()

#Semilog plot after 500 iterations
PLOT(np.arange(500,Niter+1),errors[499:],4,"Iteration Number","Log of Error",fn = pl.semilogy,title='Semilogy Plot of Error vs Iteration Number after 500 iterations',label = 'Line')
pl.legend()
pl.show()

################################################################
#Fitting a straight line to semilog plot after 500 iterations

#Ae^(Bk) for errors of large interation numbers (>= 500)
a,b = lstsq(np.c_[np.ones(Niter-499),np.arange(500,Niter+1)],np.log(errors[499:]))[0]
a = np.exp(a)
print('The values of A and B for which Ae^(Bk) fits the error after 500 iterations are:')
print(a,b)
lerr = a * np.exp(b*np.arange(500,Niter+1))

#Ae^(Bk) for entire error vector
A,B = lstsq(np.c_[np.ones(Niter),np.arange(1,Niter+1)],np.log(errors))[0]
A = np.exp(A)
print('The values of A and B for which Ae^(Bk) fits the entire error vector are:')
print(A,B)
err = A * np.exp(B*np.arange(1,Niter+1))

#Linear fitting in semilogy plot for entire error vector and for large iterations
PLOT(np.arange(1,Niter+1),errors,5,"Iteration Number","Log of Error",fn = pl.semilogy,arg3 = 'r-',title='Semilogy Plot of Error vs Iteration Number',label = 'Original Errors')
pl.semilogy(np.arange(500,Niter+1),lerr,'b-',label = 'Linearly fitted error for >500 iterations (Fit 1)')
pl.semilogy(np.arange(1,Niter+1),err,'g-',label = 'Linearly fitted error for entire error vector (Fit 2)')
pl.legend()
pl.show()

################################################################
#Generating 3-D surface plot for Potential
fig1 = pl.figure(6)
ax = p3.Axes3D(fig1)
pl.title('The 3-D Surface plot of the potential')
surf = ax.plot_surface(Y,X,Phi,rstride = 1,cstride = 1,cmap = matplotlib.cm.viridis,linewidth = 0)
fig1.colorbar(surf,shrink = 0.8,aspect = 20)
pl.show()

#Contour plot of the potential after convergence
PLOT(x,y,7,"X-axis","Y-axis",pl.contourf,Phi,"Contour Plot of Potential",cmap = matplotlib.cm.viridis)
pl.plot(X[X*X + Y*Y < 0.35*0.35],Y[X*X + Y*Y < 0.35*0.35],'ro',label = 'Points with potential 1V')
pl.legend()
pl.show()

################################################################
#Generating quiver plot of Current density

Jy = pl.zeros((Ny,Nx))
Jx = pl.zeros((Ny,Nx))
Jx[:,1:-1] = 0.5*(Phi[:,0:-2]-Phi[:,2:])
Jy[1:-1,:] = 0.5*(Phi[0:-2,:]-Phi[2:,:])
pl.figure(8)
pl.quiver(x,y,Jx,Jy)
pl.plot(X[X*X + Y*Y < 0.35*0.35],Y[X*X + Y*Y < 0.35*0.35],'ro',label = 'Points with potential 1V')
pl.xlabel('X-axis')
pl.ylabel('Y-axis')
pl.legend(loc = 1)
pl.title('Vector plot of current density')
pl.show()

################################################################
#Estimating Temperature by iterating and converging

Temp = np.zeros((Ny,Nx))*300
for k in range(Niter):
	#Interior Points
	Temp[1:-1,1:-1] = 0.25*(Temp[1:-1,0:-2]+Temp[1:-1,2:]+Temp[0:-2,1:-1]+Temp[2:,1:-1]+(Jx[1:-1,1:-1]**2)+(Jy[1:-1,1:-1]**2))
	Temp[0,:] = 300  #Boundary condition for bottom side
	Temp[1:-1,0] = Temp[1:-1,1]  #Boundary condition for left side
	Temp[1:-1,-1] = Temp[1:-1,-2]  #Boundary condition for right side
	Temp[-1,:] = Temp[-2,:]  #Boundary condition for top side
	Temp[ii] = 300  #Condition for 1V points

#Contour plot of the Temperature
PLOT(x,y,9,"X-axis","Y-axis",pl.contourf,Temp,"Contour Plot of Temperature",cmap = matplotlib.cm.hot)
pl.show()
