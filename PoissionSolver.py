import numpy as np
import time 

tic=time.time()

# number of grids in x/y-direction
Nx = 64 
Ny = 64 

# size of computational domain
Lx = 1.0
Ly = 1.0

# size of grid
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

# array for storing source
Mass           = np.zeros((Nx, Ny))

# array for storing numerical potential
SimulPotential = np.zeros((Nx, Ny))

# array for storing exact potential
ExactPotential = np.zeros((Nx, Ny))


# maximum number of iterations
MaxItr = 100000

# assign electric charge in source array
for i in range(Nx):
  for j in range(Ny):
    x = i*dx
    y = j*dy

    Mass[i][j] = 2*x*(y-1)*(y-2*x+x*y+2)*np.exp(x-y)


# Given BC condition along x-direction
for ix in range(Nx):
  SimulPotential[ix][Ny-1] = 0.0
  SimulPotential[ix][   0] = 0.0
  ExactPotential[ix][Ny-1] = 0.0
  ExactPotential[ix][   0] = 0.0


# Given BC condition along y-direction
for iy in range(Nx):
  SimulPotential[   0][iy] = 0.0
  SimulPotential[Nx-1][iy] = 0.0
  ExactPotential[   0][iy] = 0.0
  ExactPotential[Nx-1][iy] = 0.0


# initial guess for performing relaxation
#for ix in range(1, Nx-2):
#  for iy in range(1, Ny-2):
#    SimulPotential[ix][iy] = 0.5
SimulPotential[1:-2,1:-2] = 0.5




# machine epsilon in double precision
EPSILON = np.finfo(np.float64).eps

# terminate iteration when Error < Threshold
Threshold = 100*EPSILON


# number of iteration
itr = 0

#########################################
# Main Loop for relaxation
# Successive Overrelaxation method (SOR)
#########################################

while ( itr < MaxItr ):

#   the L1-norm error between the consecutive iteration cycles
    Error = 0.0

    itr += 1

    UpdatedSimulPotential = 0.25 * ( SimulPotential[1:-1,0:-2] 
                                   + SimulPotential[0:-2,1:-1]
                                   + SimulPotential[1:-1,2:  ] 
                                   + SimulPotential[2:  ,1:-1]  - dx * dy * Mass[1:-1,1:-1] )

  # sum of error
    Error = np.sum( np.absolute( np.subtract( UpdatedSimulPotential , SimulPotential[1:-1,1:-1] ) ) )
    #Error = np.sum( np.absolute(  UpdatedSimulPotential - SimulPotential[1:-1,1:-1]  ) )
  

  # calculate L-1 norm error
    Error /= ( ( Nx - 2 ) * ( Ny - 2 ) )


  # copy data from UpdatedSimulPotential to SimulPotential
    SimulPotential[1:-1,1:-1] = UpdatedSimulPotential

    if ( Error < Threshold ):
     break

# exact solution
# Ref: https://math.stackexchange.com/questions/1251117/analytic-solution-to-poisson-equation

# array for storing relative error
RelativeError = np.zeros((Nx, Ny))

# L1-norm error between exact and numerical solution
L1Error = 0.0

for i in range(1,Nx-1):
  for j in range(1,Ny-1):
     x = i*dx
     y = j*dy

     ExactPotential[i][j] = x*y*(1-x)*(1-y)*np.exp(x-y)

	 # calculate relative error between exact and numerical solution
     RelativeError[i][j] = 1 - ExactPotential[i][j] / SimulPotential[i][j]

	 # L1-norm error between exact and numerical solution
     if ( RelativeError[i][j] == RelativeError[i][j] ):
	      L1Error += abs( RelativeError[i][j] )


L1Error /= ( (Nx-2)*(Ny-2) )



# dump data to disk

f = open("potential.dat","w+")

# header
f.write("#L1Error: %20.16e\n" % L1Error)
f.write("#iterations: %d\n" % itr)
f.write("#========================================================\n")
f.write("#%19s%20s%20s%20s%20s%20s\n" %  ("x[1]", "y[2]", "Mass[3]", "Potential[4]",  "ExactPotential[5]",  "RelativeError[6]") )

#  potential data
for x in range(Nx):
  for y in range(Ny):
       f.write( "%20.7e%20.7e%20.7e%20.7e%20.7e%20.7e\n" % (  x*dx, y*dy, Mass[x][y], SimulPotential[x][y], ExactPotential[x][y], RelativeError[x][y] ) )

toc=time.time()
