import numpy as np
import time 

tic=time.time()

# number of grids in x/y-direction
Nx = 64 
Ny = 64 

# size of computational domain
Lx = 1.0
Ly = 1.0

# maximum number of iterations
MaxItr = 100000

# x coordinate of charge
DistX = 0.5*Lx
# y coordinate of charge
DistY = 0.1*Ly






# size of grid
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

# array for storing source
Charge         = np.zeros((Nx, Ny))

# array for storing numerical potential
SimulPotential = np.zeros((Nx, Ny))


# assign electric charge in source array
for i in range(Nx):
  for j in range(Ny):
    x = i*dx
    y = j*dy

    if ( abs(x - DistX ) <= dx and abs(y - DistY ) <= dy ):
      Charge[i][j] = 1.0
    else:
      Charge[i][j] = 0.0



# BC at y=0
SimulPotential[ 0:-1, 0:1] = 0.0

# BC at y=Ly
SimulPotential[ 0:-1,-1: ] = 0.0

# BC at x=0
SimulPotential[ 0: 1, 0: ] = 0.0

# BC at x=Lx
SimulPotential[-1:  , 0: ] = 0.0



# initial guess for performing relaxation
SimulPotential[1:-2,1:-2] = 0.5



# machine epsilon in double precision
EPSILON = np.finfo(np.float64).eps

# terminate iteration when Error < Threshold
Threshold = 100*EPSILON


# counter for iteration
itr = 0

#########################################
# Main Loop for relaxation
# Jacobi method
#########################################

while ( itr < MaxItr ):

#   the L1-norm error between the consecutive iteration cycles
    Error = 0.0

    itr += 1

    UpdatedSimulPotential = 0.25 * ( SimulPotential[1:-1,0:-2] 
                                   + SimulPotential[0:-2,1:-1]
                                   + SimulPotential[1:-1,2:  ] 
                                   + SimulPotential[2:  ,1:-1]  - dx * dy * Charge[1:-1,1:-1] )

  # sum of error
    Error = np.sum( np.absolute( np.subtract( UpdatedSimulPotential , SimulPotential[1:-1,1:-1] ) ) )
  

  # calculate L-1 norm error
    Error /= ( ( Nx - 2 ) * ( Ny - 2 ) )


  # copy data from UpdatedSimulPotential to SimulPotential
    SimulPotential[1:-1,1:-1] = UpdatedSimulPotential

    if ( Error < Threshold ):
     break



# dump data to disk

f = open("potential.dat","w+")

# header
f.write("#iterations: %d\n" % itr)
f.write("#========================================================\n")
f.write("#%19s%20s%20s%20s\n" %  ("x[1]", "y[2]", "Charge[3]", "Potential[4]") )

#  potential data
for x in range(Nx):
  for y in range(Ny):
       f.write( "%20.7e%20.7e%20.7e%20.7e\n" % (  x*dx, y*dy, Charge[x][y], SimulPotential[x][y] ) )

toc=time.time()
