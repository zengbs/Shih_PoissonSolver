import numpy as np


# number of grids in x/y-direction
Nx = 256
Ny = 256

# size of computational domain
Lx = 1.0
Ly = 1.0

# size of grid
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

# array for storing source
Mass = np.zeros((Nx, Ny))

# array for storing numerical potential
NumericalPotential = np.zeros((Nx, Ny))

# array for storing exact potential
ExactPotential     = np.zeros((Nx, Ny))


# assign electric charge in source array
for ix in range(Nx):
  for iy in range(Ny):
    x = ix*dx
	y = jy*dy

    Mass = 2*x*(y-1)*(y-2*x+x*y+2)*np.exp(x-y)


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
for ix in range(Nx):
  for iy in range(Ny):
    NumericalPotential[ix][iy] = 0.5




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

while ( True ):

#   the L1-norm error between the consecutive iteration cycles
    Error = 0.0

    itr++

#   update odd cells first
    for ix in range( 1, Nx-1):
      for iy in range( 1, Ny-1):
        if ( (ix+iy)%2 == 1 ):

          correction = 0.25 * w * (   SimulPotential[x+1][y  ]
                                    + SimulPotential[x-1][y  ]
                                    + SimulPotential[x  ][y+1]
                                    + SimulPotential[x  ][y-1] - dx*dy * Mass[x][y] - 4.0 * SimulPotential[x][y] )

      	# sum of error
          Error += FABS( correction/SimulPotential[x][y] )
        
      	# update grids
          SimulPotential[x][y] += correction

#   then update even cells
    for ix in range( 1, Nx-1):
      for iy in range( 1, Ny-1):
        if ( (ix+iy)%2 == 0 ):

          correction = 0.25 * w * (   SimulPotential[x+1][y  ]
                                    + SimulPotential[x-1][y  ]
                                    + SimulPotential[x  ][y+1]
                                    + SimulPotential[x  ][y-1] - dx*dy * Mass[x][y] - 4.0 * SimulPotential[x][y] )

      	# sum of error
          Error += FABS( correction/SimulPotential[x][y] )
        
      	# update grids
          SimulPotential[x][y] += correction

      	# calculate L-1 norm error
          Error /= ( Nx - 2 ) * ( Ny - 2 )

          if ( Error < Threshold ):
           break

# exact solution
# Ref: https://math.stackexchange.com/questions/1251117/analytic-solution-to-poisson-equation

# L1-norm error between exact and numerical solution
L1Error = 0.0

for ix in range(Nx):
  for iy in range(Ny):
     x = ix*dx
     y = iy*dy

     ExactPotential[i][j] = x*y*(1-x)*(1-y)*exp(x-y)

	 # calculate relative error between exact and numerical solution
     RelativeError[i][j] = 1 - ExactPotential[i][j] / Potential[i][j]

	 # L1-norm error between exact and numerical solution
	 if ( RelativeError[i][j] == RelativeError[i][j] ):
	      L1Error += FABS( RelativeError[i][j] )


L1Error /= (Nx-2)*(Ny-2)



# dump data to disk

f = open("potential.dat","w+")

# header
f.write("#L1Error: %20.16e\n" % L1Error);
f.write("#iterations: %d\n" % itr);
f.write("#========================================================\n");
f.write("%13s  %14s  %14s %16s %20s %20s\n","#x[1]", "y[2]", "Mass[3]", "Potential[4]", "ExactPotential[5]", "RelativeError[6]" );

#  potential data
for ix in range(Nx):
  for iy in range(Ny):
       f.write( "%10.7e   %10.7e   %10.7e   %10.7e   %10.7e   %10.7e\n" % (  x*dx, y*dy, Mass[x][y], Potential[x][y], ExactPotential[x][y], RelativeError[x][y] ) );
