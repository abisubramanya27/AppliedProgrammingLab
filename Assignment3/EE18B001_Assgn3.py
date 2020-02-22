####################################################
####### Applied Programming - Assignment 3 #########
####### Done by Abishek S ##########################
####################################################

#Importing libraries
import pylab as pl
import scipy.special as sp
import sys

#Part1 done in a separate Python file given

#Part2 - Extracting data from fitting.dat
try:
	data = pl.loadtxt('./fitting.dat')
#File Handling Errors
except:
	print(f'Unable to Open/Process fitting.dat')
	sys.exit()

def g(t,A = 1.05,B = -0.105):
	return A*(sp.jn(2,t)) + B*t

t = data[:,0]
Y = data[:,1:]
true_y  = g(t)

####################################################

#Part3 - Plotting the data obtained
sigma = pl.logspace(-1,-3,9)  #The Standard Deviations of noise
fig = pl.figure(0)
pl.title('Plot of Data')
pl.plot(t,pl.c_[Y,true_y])
pl.xlabel(r'$t$')
pl.legend(list(sigma) + ['true value'])
pl.grid(True)
#pl.savefig('dataplot.png')
pl.show()
pl.close(fig)

####################################################

#Part4 - Plotting the true value
fig = pl.figure(0)
pl.title('Plot of True value')
pl.plot(t,true_y,label = 'true value') #We already obtained the true value from g() in true_y
pl.xlabel(r'$t$')
pl.ylabel(r'$1.05*J(t)-0.105t$')
pl.grid(True)
pl.legend()
#pl.savefig('truevalue.png')
pl.show()
pl.close(fig)

####################################################

#Part5 - Plotting ErrorBars
fig = pl.figure(1)
pl.title('Errorbars')
pl.plot(t,pl.c_[true_y,Y[:,0]])
stdev = pl.std(Y[:,0]-true_y)
pl.errorbar(t[::5],Y[::5,0],stdev,fmt='ro')
pl.xlabel(r'$t$')
pl.legend(['True value','Noisy curve'])
#pl.savefig('errorbar.png')
pl.show()
pl.close(fig)

####################################################

#Part6 - Defining g_new based on Matrix Multiplication
def g_new(t,A = 1.05,B = -0.105):
	mat = pl.c_[sp.jn(2,t),t]
	param = pl.array([A,B])
	return pl.matmul(mat,param)

true_y_new = g_new(t)
if (true_y == true_y_new).all() == True:
	print('The two functions are equal')
else:
	print('The functions are not equal')

####################################################

#Part7 - Finding Mean Squared Error E_{i,j}
A_list = pl.linspace(0,2,21)
B_list = pl.linspace(-0.2,0,21)
#E is 2-D matrix where E[i,j] stores mean squared error for parameters A = A_list[i],B = B_list[j]
E = [[(1/len(true_y)) * sum(map(lambda x: x**2,Y[:,0]-g_new(t,a,b))) for b in B_list] for a in A_list]

####################################################

#Part8 - Plotting Contour graph
fig = pl.figure(2)
pl.title('Contour Plot')
pl.contour(A_list,B_list,E,40)
pl.plot(1.05,-0.105,'ro',label = 'Exact Value')
pl.annotate(s = 'Exact Value',xy = [0.8,-0.1])
#pl.savefig('contour.png')
pl.show()
pl.close(fig)

####################################################

#Part9 - Obtaining best estimate for A and B
(A,B),(mse),*_ = pl.linalg.lstsq(pl.c_[sp.jn(2,t),t],Y[:,0],rcond = None)
print('A : ',A)
print('B : ',B)

####################################################

#Part10 - Plotting error in estimating A and B vs noise sigma
mse_list = []
pred_A = []
pred_B = []
for i in range(Y.shape[1]):
	(a,b),(mse),*_ = pl.linalg.lstsq(pl.c_[sp.jn(2,t),t],Y[:,i],rcond = None)
	mse_list.append(mse)
	pred_A.append(a)
	pred_B.append(b)

pred_A = pl.array(pred_A)
pred_B = pl.array(pred_B)
fig = pl.figure(3)
pl.title('Error vs Noise')
#pl.plot(sigma,mse_list,label = 'Mean Squared Error')
pl.plot(sigma,pl.absolute(pred_A-1.05),'ro',label = '|Ao-Ap|')
pl.plot(sigma,pl.absolute(pred_B+0.105),'go',label = '|Bo-Bp|')
pl.legend()
pl.xlabel(r'Noise sigma')
#pl.savefig('errorvsnoise.png')
pl.show()
pl.close(fig)

####################################################

#Part11 - Plotting log of error in estimating A and B vs log of noise' sigma
fig = pl.figure(4)
pl.title('Log Error vs Log Noise')
pl.loglog(sigma,mse_list,label = 'Log of Mean Squared Error')
pl.loglog(sigma,pl.absolute(pred_A-1.05),'ro',label = 'log|Ao-Ap|')
pl.loglog(sigma,pl.absolute(pred_B+0.105),'go',label = 'log|Bo-Bp|')
pl.legend()
pl.xlabel(r'Noise sigma')
#pl.savefig('loglogplot.png')
pl.show()
pl.close(fig)

############### END OF PROGRAM ######################
