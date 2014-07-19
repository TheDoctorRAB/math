########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.19.July.2014
########################################################################
# This contains functions for matrix algebra.
# With the overall goal of computing a matrix inverse.
# This routine uses the minimum number of operations.
# This is for a sparse, square tri band matrix only.
# The matrix is positive definite by definition.
# This can be usefule for the finite element method. 
# Sometimes an inverse might be needed for FEM derivatives. 
# But can be modified without too much trouble for quad band, etc.
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
#######
#
#
#
####### efficient matrix storage
# sparse matrices contain a lot of zeroes
# these should not be stored because it wastes spaces and CPU time
# this function takes the node x node square, sparse tri band matrix
# and prepares a storage matrix of node-1 x 4 for the nonzero elements
# 00=00, 01=01, 10=02, .5(11)=03, .5(11)=10, 12=11, 21=12, etc.
###
#
###
def matrix_storage(node,main_matrix):
###
# initialize storage matrix
    storage_matrix=numpy.zeros(((node-1),4))
###
# load storage matrix
# load the 00 element
    storage_matrix[0,0]=main_matrix[0,0]
###
# load the bottom corner
    storage_matrix[node-2,1]=main_matrix[node-2,node-1]
    storage_matrix[node-2,2]=main_matrix[node-1,node-2]
    storage_matrix[node-2,3]=main_matrix[node-1,node-1]
###
# load the diagonal into the i4,i+1,0 positions by half
    for i in range(1,node-1):
        storage_matrix[i,0]=0.5*main_matrix[i,i]
        storage_matrix[i-1,3]=0.5*main_matrix[i,i]
# end i
###
# load the upper diagonal (i,1) and lower diagonal (i,2)
    for j in range(0,node-1):
        storage_matrix[j,1]=main_matrix[j,j+1]
        storage_matrix[j,2]=main_matrix[j+1,j]
# end j
###
    return(storage_matrix)
########################################################################
#
#
#
####### factorization
# the inverse requires the original matrix to be factored by LDU
# L=lower unit triangular matrix
# D=diagonal matrix
# U=upper unit triangular matrix
# however, the 1s in L and U do not need to be stored
# and only nonzeros should be stored
# the inverse of D is needed, but this is simple to compute
###
#
###
def matrix_factorization(node,storage_matrix):
###
# initialize the matrix factors (L,D,U) 
    l_factor=numpy.zeros((node-1))
    d_factor=numpy.zeros((node))
    u_factor=numpy.zeros((node-1))
# initialize inverse of D
    inverse_d=numpy.zeros((node))
###
# compute the diagonal factor
# first element
    d_factor[0]=storage_matrix[0,0]
###
# load the middle diagonals
    for i in range(1,node-1):
        d_factor[i]=(storage_matrix[i-1,3]+storage_matrix[i,0])-(storage_matrix[i-1,1]*storage_matrix[i-1,2])/d_factor[i-1]
# end i
###
# node,node element
    d_factor[node-1]=storage_matrix[node-2,3]-(storage_matrix[node-2,1]*storage_matrix[node-2,2])/d_factor[node-2]
###
# compute the upper, lower triangular
    for j in range(0,node-1):
        l_factor[j]=storage_matrix[j,2]/d_factor[j]
        u_factor[j]=storage_matrix[j,1]/d_factor[j]
###
# compute the diagonal inverse
    for k in range(0,node):
        inverse_d[k]=d_factor[k]**(-1)
# end k
###        
    return(l_factor,d_factor,u_factor,inverse_d)    
########################################################################
#
#    
#
####### inverse
# the inversion of the matrix sparse matrix is a recursive scheme
# starting with the nn element
# even though the original matrix is sparse
# the inverse will be n x n with no zero elements
###
#
###
def invert_matrix(node,l_factor,d_factor,u_factor,inverse_d):
###
# initialize inverse matrix
    inverse_matrix=numpy.zeros((node,node))
###
# compute the nn element
    inverse_matrix[node-1,node-1]=inverse_d[node-1]
###
# and the rest
    for i in range(node-2,-1,-1):
# node level row
        inverse_matrix[node-1,i]=(-1)*inverse_matrix[node-1,i+1]*l_factor[i]
        for j in range (i,-1,-1):
            if j==i:
# diagonal
                inverse_matrix[j,j]=d_factor[j]**(-1)-u_factor[j]*inverse_matrix[j+1,j]
            else:
# lower triangular by row
                inverse_matrix[i,j]=-(1)*inverse_matrix[i,j+1]*l_factor[j]
# upper triangular by column
            inverse_matrix[j,i+1]=(-1)*u_factor[j]*inverse_matrix[j+1,i+1]
# end j
# end i
###
# 
###
    return(inverse_matrix)
########################################################################
#      EOF
########################################################################