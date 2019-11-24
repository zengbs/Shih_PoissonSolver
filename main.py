import numpy as np
import time 
from PoissonSolver import *

tic=time.time()


##########################################
# input general parameters                 
##########################################
                                                 
# number of grids in x/y-direction               
Nx = 128                                          
Ny = 128
                                                 
# size of computational domain                   
Lx = 1.0                                         
Ly = 1.0                                         
                                                 
# maximum number of iterations                   
MaxItr = 100000                                  

                                                 
# size of grid
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)


##########################################
# input userdefined parameters                 
##########################################

# Problem 1 ( a charge around conductor plate )
#-----------------------------------------------------------
#f = open("charge_plate.dat","w+")
## x coordinate of charge                         
#DistX = 0.5*Lx                                   
## y coordinate of charge                         
#DistY = 0.1*Ly                                   
#                                                 
## charge
#Q = 1.0

# Problem 2 ( parallel plate capacitor )
#-----------------------------------------------------------
f = open("capacitor.dat","w+")
# distance between parallel plate
Dist = 0.2*Ly

# thickness of plate
Thickness = 2.0*dy

# length of plate
Length = 0.5*Lx

# position of center of capacitor
CenterX = 0.5*Lx
CenterY = 0.5*Ly
#
## surface charge on parallel plate
Q = 1.0

## Problem 3 ( mirror charge )
##-----------------------------------------------------------
#
#f = open("mirror_charge.dat","w+")
## x coordinate of center of charges
#DistX = 0.5*Lx                                   
#
## y coordinate of center of charges
#DistY = 0.5*Ly                                   
#
## distance betrween charges
#Dist   = 0.2*Ly                                                 
#
## charge
#Q = 1.0


##########################################
# declarations and definitions
##########################################

# array for storing source
Charge         = np.zeros((Nx, Ny))

# array for storing numerical potential
SimulPotential = np.zeros((Nx, Ny))

# array for storing exact potential
#ExactPotential = np.zeros((Nx, Ny))

# machine epsilon in double precision
EPSILON = np.finfo(np.float64).eps

# terminate iteration when Error < Threshold
Threshold = 100*EPSILON



##########################################
# assign electric charge in source array
##########################################


for i in range(Nx):
  for j in range(Ny):
    x = i*dx
    y = j*dy

# Problem 1 ( a charge around conductor plate )
#-----------------------------------------------------------
#    if ( abs(x - DistX ) <= dx and abs(y - DistY ) <= dy ):
#      Charge[i][j] = Q
#    else:
#      Charge[i][j] = 0.0


# Problem 2 ( parallel plate capacitor )
#-----------------------------------------------------------
    if   ( abs(x-CenterX) <= 0.5*Length and abs(y-CenterY-0.5*Dist) <= 0.5*Thickness ):
       Charge[i][j] = Q 
    elif ( abs(x-CenterX) <= 0.5*Length and abs(y-CenterY+0.5*Dist) <= 0.5*Thickness ):
       Charge[i][j] = -Q 
    else:
       Charge[i][j] = 0.0

# Problem 3 ( mirror charge )
##-----------------------------------------------------------
#    if   ( ( abs(x - DistX ) <= dx and abs(y - DistY - 0.5*Dist ) <= dy ) ):
#      Charge[i][j] = Q
#    elif ( ( abs(x - DistX ) <= dx and abs(y - DistY + 0.5*Dist ) <= dy ) ):
#      Charge[i][j] = -Q
#    else:
#      Charge[i][j] = 0.0
#  
#     
#
#
# Exact solution
#-----------------------------------------------------------
#    Charge[i][j] = 2*x*(y-1)*(y-2*x+x*y+2)*np.exp(x-y)


##########################################
# assign boundary conditions
##########################################

# BC at y=0
SimulPotential[ 0:-1, 0:1] = 0.0

# BC at y=Ly
SimulPotential[ 0:-1,-1: ] = 0.0

# BC at x=0
SimulPotential[ 0: 1, 0: ] = 0.0

# BC at x=Lx
SimulPotential[-1:  , 0: ] = 0.0

## BC at y=0
#ExactPotential[ 0:-1, 0:1] = 0.0
#
## BC at y=Ly
#ExactPotential[ 0:-1,-1: ] = 0.0
#
## BC at x=0
#ExactPotential[ 0: 1, 0: ] = 0.0
#
## BC at x=Lx
#ExactPotential[-1:  , 0: ] = 0.0


##########################################
# initial guess for performing relaxation
##########################################

SimulPotential[1:-2,1:-2] = 0.5


##########################################
# perform relaxation
##########################################

SimulPotential, itr = PoissionSolver( SimulPotential, Charge, Nx, dx, Ny, dy, Threshold, MaxItr )


##########################################
# compare numerical solution with exact solution
##########################################
## exact solution
## Ref: https://math.stackexchange.com/questions/1251117/analytic-solution-to-poisson-equation
#
## array for storing relative error
#RelativeError = np.zeros((Nx, Ny))
#
## L1-norm error between exact and numerical solution
#L1Error = 0.0
#
#for i in range(1,Nx-1):
#  for j in range(1,Ny-1):
#     x = i*dx
#     y = j*dy
#
#     ExactPotential[i][j] = x*y*(1-x)*(1-y)*np.exp(x-y)
#
#	 # calculate relative error between exact and numerical solution
#     RelativeError[i][j] = 1 - ExactPotential[i][j] / SimulPotential[i][j]
#
#	 # L1-norm error between exact and numerical solution
#     if ( RelativeError[i][j] == RelativeError[i][j] ):
#	      L1Error += abs( RelativeError[i][j] )
#
#
#L1Error /= ( (Nx-2)*(Ny-2) )

##########################################
# dump data to disk
##########################################


# header
f.write("#iterations: %d\n" % itr)
f.write("#========================================================\n")
f.write("#%19s%20s%20s%20s\n" %  ("x[1]", "y[2]", "Charge[3]", "Potential[4]") )

# potential data
for x in range(Nx):
  for y in range(Ny):
       f.write( "%20.7e%20.7e%20.7e%20.7e\n" % (  x*dx, y*dy, Charge[x][y], SimulPotential[x][y] ) )




##########################################
# dump exact soultion data to disk
##########################################

## header
#f.write("#L1Error: %20.16e\n" % L1Error)
#f.write("#iterations: %d\n" % itr)
#f.write("#========================================================\n")
#f.write("#%19s%20s%20s%20s%20s%20s\n" %  ("x[1]", "y[2]", "Charge[3]", "Potential[4]",  "ExactPotential[5]",  "RelativeError[6]") )
#
##  potential data
#for x in range(Nx):
#  for y in range(Ny):
#       f.write( "%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e\n" % (  x*dx, y*dy, Charge[x][y], SimulPotential[x][y], ExactPotential[x][y], RelativeError[x][y] ) )


toc=time.time()
