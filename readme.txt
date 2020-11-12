This is a series of program written in Fortran to reproduce the result of Tully's Fewest Switches Surface Hopping model 
and some variant (phase-corrected surface hopping and branching-corrected surface hopping) to make the results closer to 
the full quantum one. 

For a specifically given initial momentum, the program will produce the final probability of each reaction channel 
(transmission to the upper surface, transmission to the lower surface or reflection to the lower surface). The simulation
 time step dt should be chosen sufficiently small to make sure the results are already converge, and the number of trajectories
 is also need to be large enough to remove the statistical error.