import numpy as np
import os.path
import matplotlib.pyplot as plt
import scipy.sparse as sparse

DIM = 1
p = 6
res = 2**p
degree_max = 10 #number of basis along one dimension

def mult_by_monomial( alpha ):
    #produces multiplication operator for monomial
    global DIM, degree_max
    out = sparse.eye( degree_max, k=alpha[0], format='dia')
    for d in range(1,DIM):
        store = sparse.eye( degree_max, k=alpha[d], format='dia')
        out = sparse.kron( store , out )
    return out

def basis_function(beta,x):
    #one dimesional basis function
    X = (np.abs(x) < 7)*x
    if beta > 0:
        return ( X*np.exp( - X**2 / (2.*beta) ) )**beta
    else:
        return np.exp( - X**2 / 2.0 )

def derivative_of_basis_function(beta,x):
    if beta > 0:
        return beta * psi(beta-1,x) - psi(beta+1,x)
    else:
        return -psi(beta+1,x)

def translation_gen(d):
    # generates a one-dimensional translation
    global degree_max,DIM
    data = -np.ones( [2 , degree_max] )
    data[0,:] = np.arange(0, degree_max)
    offsets = np.array( [ 1 , -1 ] , dtype = int)
    ide = sparse.eye(degree_max,format='dia')
    ddx = sparse.dia_matrix( (data,offsets) , shape = (degree_max,degree_max) )
    if d==0:
        out = ddx
    else:
        out = ide 
    for i in range(1,DIM):
        if d==i:
            out = sparse.kron(out,ddx)
        else:
            out = sparse.kron(out,ide)
    return out

x_arr = np.linspace(-3,3,res)
X,Y = np.meshgrid(x_arr,x_arr) 

#fname = "basis_array_DIM_{:d}.npy".format(DIM)
#if os.path.isfile(fname):
#    basis_array = np.load(fname)
#    print "Loaded basis array"
#else:
#    if DIM ==1:
#        basis_array = np.zeros([degree_max,res])
#    elif DIM==2:
#        basis_array = np.zeros([degree_max**2,res,res])
#    elif DIM==3:
#        basis_array = np.zeros([degree_max**3,res,res,res])
#    for beta in range(0,10):
#        basis_array[beta] = basis_function( beta , X)
#    np.save( fname , basis_array)
#    print "Created basis array"

#DESCRIBE YOUR VECTOR_FIELD HERE
f = [{} , {}]
f[0][(1,0)] = -1.0
f[1][(0,0)] = 1.0

# This is the vector field dx/dt = -x , dy/dt = 1
Operator = translation_gen(1)
psi_initial = np.zeros( degree_max )
psi_initial[0] = 1.0
psi_initial[1] = 0

monomial_indices  = [ (a_1,a_2) for a_1 in range(0,degree_max)
            for a_2 in range(0,degree_max) ]

print monomial_indices


