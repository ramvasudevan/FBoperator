import numpy as np
import os.path
import matplotlib.pyplot as plt
import scipy.sparse as sparse

DIM = 1
p = 6
res = 2**p
n_basis = 10 #number of basis along one dimension

def mult_by_monomial( alpha ):
    #produces multiplication operator for monomial
    global DIM, n_basis
    out = sparse.eye( n_basis, k=alpha[0], format='dia')
    for d in range(1,DIM):
        store = sparse.eye( n_basis, k=alpha[d], format='dia')
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
    global n_basis,DIM
    data = -np.ones( [2 , n_basis] )
    data[0,:] = np.arange(0, n_basis)
    offsets = np.array( [ 1 , -1 ] , dtype = int)
    ide = sparse.eye(n_basis,format='dia')
    ddx = sparse.dia_matrix( (data,offsets) , shape = (n_basis,n_basis) )
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


X = np.linspace(-3,3,res)

if os.path.isfile("basis_array.npy"):
    basis_array = np.load("basis_array.npy")
    print "Loaded basis array"
else:
    basis_array = np.zeros([n_basis,res])
    for beta in range(0,10):
        basis_array[beta] = basis_function( beta , X)
    np.save( "basis_array" , basis_array)
    print "Created basis array"

#DESCRIBE YOUR VECTOR_FIELD HERE
f = [{} , {}]
f[0][(1,0)] = -1.0
f[1][(0,0)] = 1.0

# This is the vector field dx/dt = -x , dy/dt = 1

Operator = translation_gen()
mult_by_x = multiply_by_x(2)
psi_initial = np.zeros( n_basis )
psi_initial[0] = 1.0
psi_initial[1] = 0

plt.plot(X , np.dot(psi_initial,basis_array) )
plt.plot(X , np.dot( mult_by_x.dot(psi_initial)  , basis_array) )
plt.grid(True)
plt.show()

