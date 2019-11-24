import numpy as np

def PoissionSolver ( SimulPotential, Charge, Nx, dx, Ny, dy, Threshold, MaxItr ):

   # counter for iteration
   itr = 0

   while ( itr < MaxItr ):
   
       # the L1-norm error between the consecutive iteration cycles
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

   return SimulPotential, itr
