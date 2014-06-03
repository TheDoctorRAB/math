########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.03.June.2014
########################################################################
# Analytical solution to the diffusion reaction equation.
# Cartesian coordinates are used. 
# This is for a semi-infinite domain with Dirichlet boundary conditions.
# The initial condition is zero.
# The solution is normalized to the initial condition.
# This is also used as a comparison for the companion FEM solution. 
########################################################################
#
########################################################################
# e,i,j,k,n are used for loop indices
########################################################################
#
#
#
####### imports
import numpy
import scipy.special
import side_conditions
import fem_functions as fem_f
import physical_constants as phyc
import make_plot
#######
#
####### read in the input data
# boundary conditions and time
bc1,bc2,time=side_conditions.side_conditions()
###
# physical constants
decay_constant,diffusion_coefficient,porosity,retardation=phyc.physical_constants()
#######
#
####### construct spatial domain
# this uses a function created for the FEM solution
# h and element are not needed for this routine
# they are needed for the FEM solution and this function just does all that
#
spatial_mesh,h,node,element=fem_f.spatial_mesh_discretization()
#######
#
####### initialize the solution matrix
# the spatial mesh is arbitrarily chosen to be 1001 points
solution=numpy.zeros((node,2))
#######
#
####### load in the spatial domain to the solution matrix
for i in range(0,node):
    solution[i,0]=spatial_mesh[i]
# end i
#######
#
#######
# compute solution
# the arguments in the exp and erfc are complicated
# so just simplify them here for coding ease
###
#
###
erfc1=numpy.sqrt((retardation)/(4*porosity*diffusion_coefficient*time))
erfc2=decay_constant*time
exp=numpy.sqrt((retardation*decay_constant)/(porosity*diffusion_coefficient))
###
# assign boundary conditions
solution[0,1]=bc1
solution[node-1,1]=bc2
for j in range(1,node-1):
    solution[j,1]=(0.5)*((scipy.exp(solution[j,0]*exp)*scipy.special.erfc((solution[j,0]*erfc1)+(erfc2))+(scipy.exp(-solution[j,0]*exp)*scipy.special.erfc((solution[j,0]*erfc1)-(erfc2)))))
    print solution[j,1],j
# end j
#######
#
####### plot
make_plot.makeplot(solution,1)
#######
# save solution
numpy.savetxt('solution_diffusion.reaction.equation.out',solution)
#######
#
#
#
########################################################################
#      EOF
########################################################################
