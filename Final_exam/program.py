########################################################
####### Applied Programming - Final Assignment #########
####### Done by Abishek S ##############################
########################################################

# Importing libraries
import pylab as pl
import argparse
from scipy.linalg import lstsq
import sys
import matplotlib

################################################################
# Getting optional arguments and setting certain parameter values for the problem
''' 
    We already know the following information for the problem :
        Lx : 10 cm - Physical length of tank along x direction
        Ly : 20 cm - Physical length of tank along y direction
        e_r : 2 - Relative permttivity of fluid in tank
'''
Lx,Ly,e_r = 0.1, 0.2, 2.0  # Lx, Ly in metres

# Using argparse to get the following parameters as optional command line arguments
parser = argparse.ArgumentParser(description='Solving for potential of a tank in the given setup  ***** Make sure distance between nodes along x and along y are equal while choosing M and N *****')
# M : The number of nodes along y, including the boundary nodes 
parser.add_argument('--M',dest='M',type=int,default=41,help = 'The number of nodes along y, including the boundary nodes (>= 2)')
# N : The number of nodes along x, including the boundary nodes
parser.add_argument('--N',dest='N',type=int,default=21,help = 'The number of nodes along x, including the boundary nodes) (>= 2)')
# delta : The desired accuracy for the potential obtained
parser.add_argument('--delta',dest='delta',type=float,default=1e-8,help = 'The desired accuracy for the potential obtained')
# Ni : The maximum number of iterations to complete for convergence
parser.add_argument('--Ni',dest='NIter_max',type=int,default=3000,help = 'The maximum number of iterations to complete for convergence')

args = parser.parse_args()
M,N,delta,NIter_max = args.M,args.N,args.delta,args.NIter_max

# Checking if distance between a node is same along x and along y
if (Ly / (M-1)) != (Lx / (N-1)):
    print('\nMake sure ditance between nodes along x and along y are equal for given M and N')
    sys.exit()

# Distance between nodes (same along x and along y) in metres
dist = Ly / (M-1)

################################################################
# Utility function for plotting data

def PLOT(x,y,label_x = r'X$\rightarrow$',label_y = r'Y$\rightarrow$',fn = pl.plot,arg3 = '-',title = "Plot",fig_no = 0,grids = True,label = '',cmap = matplotlib.cm.jet):
    '''
        Utility function for making the more repeated plots
        Takes in -
            x : Data points for x axis
            y  : Data points for y axis
            label_x : Label for x axis
            label_y : Label for y axis
            fn : Which plot function to use
            arg3 : 3rd  argument to the function - (the matrix for contour plot, the line style for normal plot)
            title : Title for the plot
            fig_no : Figure number for the plot
            grids : True is grids need to be present on the plot, False otherwise
            label : Legend label for the plot drawn
            cmap : Colour map to use for the contour plot
    '''
    pl.figure(fig_no)
    if fn == pl.contourf:
        fn(x,y,arg3,cmap = cmap)
        pl.colorbar()
    else:
        if label == '':
            fn(x,y,arg3)
        else:
            fn(x,y,arg3,label = label)
    pl.xlabel(label_x,size=15)
    pl.ylabel(label_y,size=15)
    pl.title(title,size=16)
    pl.grid(grids)

################################################################
# Function to solve Laplcae equation for given set of parameters - Part (d)

def solveLaplace(M,N,dist,k,delta,Ni_max):
    ''' 
        The function solves Laplace equation given -
            M : The number of nodes along y, including the boundary nodes
            N : The number of nodes along x, including the boundary nodes
            dist : Distance between nodes (same along x and along y)
            k : The height given as the index k corresponding to h
            delta : The desired accuracy for the potential obtained
            Ni_max : The maximum number of iterations to complete for convergence
        -----------------------------------------------
        The function returns -
            phi[M,N] : The array of solved potential values correct to delta
            Ni : Number of iterations actually carried out
            err[Ni] : The vector of errors
    '''
    # phi is a matrix where (0,0) corresponds to the lower left corner of tank and (M-1,N-1) corresponds to the top right corner
    phi = pl.zeros(shape = (M,N)) # Initialising Potential grid to zero at all nodes
    phi[-1,1:-1] = 1 # Top boundary points are at 1V
    # The bottom and side boundary points are at 0 V (grounded) which is already satisfied
    errors = pl.zeros(Ni_max)
    Ni = 0

    for iter_no in range(Ni_max):
        Ni += 1
        oldPhi = phi.copy()
        ##### Using vectorised code to execute code faster #####
        # Interior Points from 0 to (k-1)th index - average of surrounding points
        phi[1:k,1:-1] = 0.25 * (oldPhi[0:k-1,1:-1] + oldPhi[2:k+1,1:-1] + oldPhi[1:k,0:-2] + oldPhi[1:k,2:])
        # Interior Points from (k+1)th to Nth index - average of surrounding points
        phi[k+1:-1,1:-1] = 0.25 * (oldPhi[k:-2,1:-1] + oldPhi[k+2:,1:-1] + oldPhi[k+1:-1,0:-2] + oldPhi[k+1:-1,2:])
        # Interior Points at kth index - slightly different updation to handle Dn continuity
        global e_r
        phi[k,1:-1] = (e_r*oldPhi[k-1,1:-1] + oldPhi[k+1,1:-1]) / (1 + e_r)

        # The top, bottom and side boundaries are at constant potentials as initialised and unaffected by above update
        # Hence not running below code and increase the time taken
        # phi[0,:] = 0  # Bottom side
        # phi[:,0] = 0  # Left side
        # phi[:,-1] = 0  #Right side
        # phi[-1,1:-1] = 1  # Top side
        
        errors[iter_no] = (abs(phi-oldPhi)).max()
        if iter_no > 500 and errors[iter_no] < delta:
            # Running minimum 500 iterations to get a fitting for the error after 500 iterations
            # Exiting the iterations on reaching desired accuracy (i.e) when error goes below the required limit
            break
    
    errors = errors[:Ni]
    return phi,Ni,errors

################################################################
# Running the laplace equation solver for differernt h values - part (e)

EPS = 1e-9 # A constant to compare equal decimals which would differ slightly due to machine precision
e_o = 8.854e-12 # Absolute Pe_rmittivity of free space - constant


def ElectricField(phi):
    '''
        Function to calculate electric field along x and y directions using one sided derivative
        Takes in -
            phi : The potential grid
        -----------------------------------------------
        Returns - 
            Ex : Electric field along x direction {Ex[m,n] for points (m,n-0.5) - since one sided derivative used}
            Ey : Electric field along y direction {Ey[m,n] for points (m-0.5,n) - since one sided derivative used}
    '''
    global dist
    Ex = pl.zeros((M,N))
    Ey = pl.zeros((M,N))

    # Electric field is calculated as numerical one-sided derivative of potential along respective directions
    # Ex[i,j+1] = -(phi[i,j+1] - phi[i,j]) / dist 
    # and Ey[i+1,j] = -(phi[i+1,j] - phi[i,j]) / dist
    # where dist = delta x = delta y
    Ex[:,1:] = -(phi[:,1:] - phi[:,:-1]) / dist
    Ey[1:,:] = -(phi[1:,:] - phi[:-1,:]) / dist

    return Ex,Ey

def Charge_Top_Fluid(Ex,Ey):
    '''
        Function to caclulate Qtop and Qfluid
        Takes in -
            Ex : Electric field along x direction at (m,n+0.5)
            Ey : Electric field along y direction at (m+0.5,n)
        -----------------------------------------------
        Returns -
            Qtop : Charge at the top surface
            Qfluid : Charge at the surfaces in contact with the fluid
    '''
    # En_top is normal electric field (along +y direction - outward normal) at the top surface (i.e) (M-0.5,n)
    En_top = Ey[-1,:]
    # En_lside is normal electric field (along -x direction - outward normal) at the left side till fluid is present (i.e) (m,0.5)
    En_lside = -Ex[:k+1,1]
    # En_rside is normal electric field (along +x direction - outward normal) at the right side till fluid is present (i.e) (m,N-0.5)
    En_rside = Ex[:k+1,-1]
    # En_bottom is normal electric field (along -y direction - outward normal) at the bottom surface (i.e) (0.5,n)
    En_bottom = -Ey[1,:]

    global e_o,e_r,dist

    # Qtop consists of only the top wall
    Qtop = -e_o * sum(En_top) * dist  # dist is constant over summation hence brought out 

    # Qfluid consists of the side walls till height h and the bottom wall
    Qfluid = -e_o*e_r * (sum(En_lside) + sum(En_rside) + sum(En_bottom)) * dist  # dist is constant over summation hence brought out 

    return Qtop,Qfluid


# x and y are the axes of the grid to make the contour plot
x = pl.linspace(0,Lx,N)*100  # in cm
y = pl.linspace(0,Ly,M)*100  # in cm

# Qtop is the charge on the top plate
Qtop = pl.zeros(9)
# Qfluid is the charge on the walls of tank in contact with the fluid
Qfluid = pl.zeros(9)

for ind,hbyLy in enumerate([x*0.1 for x in range(1,10)]):
    # h/Ly = k/(M-1) but k and M are integers => (h*(M-1)/Ly) should be an integer
    k = (hbyLy*(M-1))
    if abs(k - int(k)) > EPS:
        # k is not an integer
        print("\nFor the given value of M and h, there doesn't exist an index corresponding to the fluid top boundary")
        sys.exit()
    k = int(k)
    phi,Ni,errors = solveLaplace(M,N,dist,k,delta,NIter_max)
    # Fitting an exponential to the error obtained during each iteration
    # Obtaining A and B by using lstsq on log(y) = log(A) + B.x , where y is error and x is the iteration number
    A,B = lstsq(pl.c_[pl.ones(Ni-499),pl.arange(500,Ni+1)],pl.log(errors[499:]))[0]
    A = pl.exp(A)
    print('\nThe values of A and B for which Ae^(Bk) fits the iteration error vector (for h/Ly = {:.1f}):'.format(hbyLy))
    print(A,B)
    print('The maximum error on extrapolating the error to infinity (for h/Ly = {:.1f}):'.format(hbyLy))
    print(-A/B * pl.exp(B*(Ni+0.5)))
    # The exponential which fits the iteration error
    error_fit = A * pl.exp(B*pl.arange(500,Ni+1))

    # Uncomment the below lines to view the potential contour plots and semilog plot of error with the extrapolation
    # PLOT(pl.arange(1,Ni+1),errors,r'$Number\ of\ iterations\rightarrow$',r'$Error\rightarrow$',pl.semilogy,'b-',r'$Log\ (error)\ vs\ Iteration\ number\ for\ h/L_y\ =\ {:.1f}$'.format(hbyLy),0,True,'True iteration error')
    # pl.semilogy(pl.arange(500,Ni+1),error_fit,'g-',label = 'Fitted error')
    # pl.legend()
    # pl.show()
    # PLOT(x,y,r'$X\ (in\ cm)\ \rightarrow$',r'$Y\ (in\ cm)\ \rightarrow$',pl.contourf,phi,r'$Contour\ plot\ of\ Potential\ for\ h/L_y\ =\ {:.1f}$'.format(hbyLy),1,cmap = matplotlib.cm.plasma)
    # pl.show()

    '''
        Since the walls of the tank are all made of conductors, we can use the fact that -D.n = sigma
        where n - the outward normal to the wall at a point and sigma - the charge density (charge/unit area) at that point

        Hence we can sum the sigma multiplied by dist (numerically equal to the small length over which sigma is the charge density) at each point to get charge per unit depth of the tank
        Q = charge over unit depth (depth of tank (Lz) is unknown and constant so we can assume it is unity(1 metre) )
    '''
    # D.n can be calculated as -e dV/dn where e is absolute permittivity of medium = e_o * e_r, dV/dn is derivative of potential along outward normal direction
    # We use numerical approximation on derivative => dv = V[i+1] - V[i] and dn = dist = delta x = delta y

    Ex,Ey = ElectricField(phi)
    Qtop[ind],Qfluid[ind] = Charge_Top_Fluid(Ex,Ey)

################################################################
# Plotting Qtop and Qfluid vs h plots - part (e)

PLOT(pl.arange(0.1,1,0.1)*Ly*100,Qtop,r'$h\ (in\ cm)\ \rightarrow$',r'$Q_{top}\ (in\ C)\ \rightarrow$',pl.plot,'-o',r'$Q_{top}\ Vs\ h$',2,True,r'$Q_{top}$')
pl.legend()
pl.show()
PLOT(pl.arange(0.1,1,0.1)*Ly*100,Qfluid,r'$h\ (in\ cm)\ \rightarrow$',r'$Q_{fluid}\ (in\ C)\ \rightarrow$',pl.plot,'-o',r'$Q_{fluid}\ Vs\ h$',3,True,r'$Q_{fluid}$')
pl.legend()
pl.show()
# Tried semilogy plot to check if Qtop vs h follows a linear trend in semilogy
# PLOT(pl.arange(0.1,1,0.1)*Ly*100,Qtop,r'$h\ (in\ cm)\ \rightarrow$',r'$log(Q_{top})\ (Q_{top}\ in\ C)\ \rightarrow$',pl.semilogy,'-o',r'$log(Q_{top})\ Vs\ h$',4,True,r'$log(Q_{top})$')
# pl.show()
# Tried a simple hyperbolic plot to check if 1/Qtop varies linearly with h
# PLOT(pl.arange(0.1,1,0.1)*Ly*100,1/(Qtop),r'$h\ (in\ cm)\ \rightarrow$',r'$1/Q_{top}\ (Q_{top}\ in\ C)\ \rightarrow$',pl.plot,'-o',r'$1/Q_{top}\ Vs\ h$',5,True,r'$1/Q_{top}$')
# pl.show()

################################################################
# Finding Ex and Ey at centre of mesh cells and showing Continuity of Dn at m = k - part (f)

def ElectricField_Centre(phi):
    '''
        Function to calculate Ex and Ey at centre of mesh cells (i.e) at (m+0.5,n+0.5)
        Takes in -
            phi : The potential grid
        -----------------------------------------------
        Returns -
            Ex_centre : Ex at centre of mesh cells
            Ey_centre : Ey at centre of mesh cells
    '''
    # Getting Ex and Ey
    Ex,Ey = ElectricField(phi)

    # Ex[m,n+1] finds Ex at (m,n+0.5). Ex at (m+0.5,n+0.5) would be -(phi @ (m+0.5,n+1) - phi @ (m+0.5,n)) / dist
    # Since we don't know phi @ (m+0.5,n), we can approximate that as average of phi @ (m,n) and phi @ (m+1,n)
    # Similarly for phi @ (m+0.5,n+1)
    # On rearranging the terms, we get Ex_centre[m+1,n+1] as average of Ex[m+1,n+1] and Ex[m,n+1], which will be Ex at (m+0.5,n+0.5)
    Ex_centre = pl.zeros((M,N))
    Ex_centre[1:,1:] = 0.5 * (Ex[1:,1:] + Ex[:-1,1:])

    # Similarly Ey_centre[m+1,n+1] is average of Ey[m+1,n+1] and Ey[m+1,n], which will be Ey at (m+0.5,n+0.5)
    Ey_centre = pl.zeros((M,N))
    Ey_centre[1:,1:] = 0.5 * (Ey[1:,1:] + Ey[1:,:-1])

    return Ex_centre,Ey_centre


# We need k for h/Ly = 0.5
k = int(0.5 * (M-1))  # k is centre along y axis
phi,Ni,errors = solveLaplace(M,N,dist,k,delta,NIter_max)

# Uncomment the below part to look at contour plot of potential for h/Ly = 0.5 case
# PLOT(x,y,r'$X\ (in\ cm)\ \rightarrow$',r'$Y\ (in\ cm)\ \rightarrow$',pl.contourf,phi,r'$Contour\ plot\ of\ Potential\ for\ h/L_y\ =\ 0.5$',6,cmap = matplotlib.cm.plasma)
# pl.show()

# Uncomment the below part to look at quiver plot of electric fields in the interior points of the grid calculated using double sided derivative
# Ex_d = pl.zeros((M,N))
# Ey_d = pl.zeros((M,N))
# Ex_d[1:-1,1:-1] = -(phi[1:-1,2:] - phi[1:-1,:-2]) / dist
# Ey_d[1:-1,1:-1] = -(phi[2:,1:-1] - phi[:-2,1:-1]) / dist
# pl.figure(7)
# pl.quiver(x,y,Ex_d,Ey_d)
# pl.xlabel(r'$X\ (in\ cm)\ \rightarrow$',size=15)
# pl.ylabel(r'$Y\ (in\ cm)\ \rightarrow$',size=15)
# pl.title('Electric field at interior points of grid',size=16)
# pl.show()


####### Ex, Ey at (m+0.5,n+0.5) - at centre of mesh cells
Ex_centre,Ey_centre = ElectricField_Centre(phi)

# Uncomment the below part to print the Ex and Ey at (m+0.5,n+0.5)
# print('\n\nEx values at (m+0.5,n+0.5) :')
# print(Ex_centre[1:,1:])

# print('\nEy values at (m+0.5,n+0.5) :')
# print(Ey_centre[1:,1:])

# x and y are the axes of the grid to make the quiver plot - coordinates of midpoint of the mesh cells
x_c = pl.arange(dist/2,Lx,dist)*100  # in cm
y_c = pl.arange(dist/2,Ly,dist)*100  # in cm
pl.figure(8)
pl.quiver(x_c,y_c,Ex_centre[1:,1:],Ey[1:,1:])
pl.xlabel(r'$X\ (in\ cm)\ \rightarrow$',size=15)
pl.ylabel(r'$Y\ (in\ cm)\ \rightarrow$',size=15)
pl.title('Electric field at centre of mesh cells',size=16)
pl.show()


# Checking continuity of Dn at m = k
# Dn = e * En, where e = absolute permittivity of medium and En = Ey at the fluid surface m = k
# We are basically checking Dn just above and below are same => (Left hand limit = Right hand limit) => continuity
# I haven't multiplied by e_o since we are finding relative percentage difference which is a ratio

# Ey just above m = k (i,e) (k+0.5,n) => Ey_centre[k+1,1:]; Ey just below m = k (i,e) (k-0.5,n) => Ey_centre[k,1:]
per_diff = abs((Ey_centre[k+1,1:]-e_r*Ey_centre[k,1:]) / Ey_centre[k+1,1:]).max() * 100  # percentage difference
print(f'\nThe relative difference between Dn just above and below is nearly : {per_diff} %')
if per_diff < 0.1:
    print('Dn is continuous at m = k')
else :
    print('Dn is not continuous at m = k')

################################################################
# Change in angle of electric field at m = k - part (g)

# Using Ex_centre and Ey_centre from above part which calculated Electric fields along x and y axes at centre of mesh cells
# Angle made by electric field with normal (y-axis) at points just below m = k (i.e) (k-0.5,n)
angle_b = pl.arctan(pl.divide(Ex_centre[k,1:],Ey_centre[k,1:]))  # Theta1
# Angle made by electric field with normal (y-axis) at points just above m = k (i.e) (k+0.5,n)
angle_t = pl.arctan(pl.divide(Ex_centre[k+1,1:],Ey_centre[k+1,1:]))  # Theta2

# Plot of change in angle of electric field (i.e) angle_t - angle_b
PLOT(x_c,angle_t-angle_b,r'$X\ (in\ cm)\ \rightarrow$',r'$Change\ in\ angle\ (in\ radians)\ \rightarrow$',pl.plot,'-o',"Change in angle of Electric Field at m = k",9)
pl.show()
# Uncomment the below line to print the change in angle of Electric field at m = k
# print(angle_t-angle_b)

# Leaving the constants behind since we are considering the ratio only
ratio =  (e_r**0.5 * pl.sin(angle_b)) / pl.sin(angle_t)  # n1.sin(Theta1) / n2.sin(Theta2)
PLOT(x_c,ratio,r'$X\ (in\ cm)\ \rightarrow$',r'$\frac{n_1sin(\theta_1)}{n_2sin(\theta_2)}\ \rightarrow$',pl.plot,'-o',"Snell's law validity",10)
pl.show()
mean_ratio = ratio.mean()
std_ratio = ratio.std()
print(f'\nThe value of n1*sin(Theta1) / n2*sin(Theta2) has mean : {mean_ratio} and standard deviation : {std_ratio}')
# The ratio has to be closer to 1 at all points for snell's law to be valid, hence standard deviation of the ratio too has to be small apart from mean of the ratio being near to 1
if abs(mean_ratio-1.0) < 0.1 and std_ratio < 0.2:
    print("Snell's law is followed")
else :
    print("Snell's law is not followed")


############################## END OF PROGRAM ##################################