import numpy as np
from scipy.optimize import linprog
from matplotlib import pyplot as plt

Nx = 20
Nt = 20
MIN_X = -1.0
MAX_X = 1.0
MIN_T = 0.0
MAX_T = 1.0
span_x = (MIN_X,MAX_X,Nx)
span_t = (MIN_T,MAX_T,Nt)
x = np.linspace( *span_x )
t = np.linspace( *span_t )
dx = x[1]-x[0]
dt = t[1]-t[0]

def derivative_op(x_min,x_max,N):
    #produces a one-dimesional derivative operator
    dx = (x_max - x_min)/float(N-1)
    out = (np.eye(N)-np.eye(N,k=-1))/dx
    out[0,1] = 1./dx
    out[0,0] = -1./dx
    out[N-1,N-1] = 1./dx
    out[N-1,N-2] = -1./dx
    return out

def unit_test_derivative_op(min_x,max_x,N):
    x = np.linspace(min_x,max_x,N)
    y = np.sin(x)
    ddx = derivative_op(min_x,max_x,N)
    dy_dx = np.dot(ddx,y)
    plt.plot( x , dy_dx-np.cos(x) , 'r-')
    plt.title('error')
    plt.show()
    return 0

#unit_test_derivative_op( *span_x)
#quit()

from scipy.linalg import kron

partial_x = kron( np.eye(Nt), derivative_op( *span_x ) )
partial_t = kron( derivative_op( *span_t ), np.eye(Nx) )

def unit_test_partial_t():
    x = np.linspace( *span_x )
    t = np.linspace( *span_t )
    #consider the function x*(t**2)
    f = np.kron( t**2 , x )
    df_dt_computed = np.dot( partial_t , f )
    df_dt_computed.resize( Nt,Nx)
    df_dt = np.kron( 2*t , x )
    df_dt.resize(Nt,Nx)
    plt.subplot(2,1,1)
    plt.imshow( df_dt )
    plt.ylabel('x')
    plt.subplot(2,1,2)
    plt.imshow( df_dt_computed )
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()
    return 0

def unit_test_partial_x():
    x = np.linspace( *span_x )
    t = np.linspace( *span_t )
    #consider the function exp(t)*sin(x)
    f = np.kron( np.exp(t) , np.sin(2*np.pi*x) )
    df_dx_computed = np.dot( partial_x , f )
    df_dx_computed.resize( Nt,Nx)
    df_dx = np.kron( np.exp(t) , 2*np.pi*np.cos(2*np.pi*x) )
    df_dx.resize(Nt,Nx)
    plt.subplot(2,1,1)
    plt.imshow( df_dx )
    plt.ylabel('x')
    plt.subplot(2,1,2)
    plt.imshow( df_dx_computed )
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()
    return 0


k = 1
from numpy.linalg import matrix_power
grad_k = matrix_power( partial_x , k )

L = partial_t + 0.0*partial_x

#Create restriction operators
#X_T = [-.5 , 0.5]
XT_indices = [ i for i in xrange(Nx) if (i > Nx/4) and (i < 3*Nx/4) ]
XTc_indices = [i for i in xrange(Nx) if i not in XT_indices ]

Ix = np.eye(Nx)
at_time_T = np.zeros( Nt )
at_time_T[Nt-1] = 1.0
restrict_XT = np.kron(at_time_T, Ix[ XT_indices,:])
restrict_XTc = np.kron(at_time_T,  Ix[XTc_indices,:] )

on_left = np.zeros(Nx)
on_left[0] = 1.0
boundary_condition = np.kron( np.eye(Nt) , on_left )

#linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None, method='simplex', callback=None, options=None)
def unit_test_restriction():
    rho = np.random.randn( Nt, Nx)
    rho_T = rho[Nt-1,:]
    rho_T_restricted = np.dot(restrict_XT , rho.flatten())
    import matplotlib.pyplot as plt
    plt.plot( x , rho_T, 'b')
    plt.plot( x[XT_indices] , rho_T_restricted , 'ro')
    plt.show()
    return (rho_T[XT_indices]==rho_T_restricted).all()

#cost function is the integral at time 0
delta_0 = np.zeros(Nt)
delta_0[0] = 1.0

#  All operators act on (w(x) ,v(t,x))

#c(w,v) = int_X w(x) dx
c = np.zeros( Nx + Nx*Nt )
c[0:Nx] = np.ones(Nx)*dx 

A_list = []
b_list = []

#Louiville constraint
A_list.append( np.hstack( [np.zeros((L.shape[0],Nx)), L ]) )
b_list.append( np.zeros(L.shape[0]) )

# w(x) >= v(0,x) + 1, or v(0,x) - w(x) <= -1
evaluate_at_0 = np.kron( delta_0 , np.eye(Nx) )
A_list.append( np.hstack( [-np.eye(Nx) , evaluate_at_0] ) )
b_list.append( -np.ones(Nx) )

#v(T,x) >=0 for x in X_T, or -v(T,x) <= 0 for x in X_T
A_list.append( np.hstack( [np.zeros( (restrict_XT.shape[0],Nx) ) , -restrict_XT] ) )
b_list.append( np.zeros( restrict_XT.shape[0] ) )

print restrict_XT.shape

#w(x) >=0 for all x in X or -w(x) <= 0
A_list.append( np.hstack([ -np.eye(Nx) , np.zeros( (Nx,Nx*Nt) )] ) )
b_list.append( np.zeros( Nx ) )

A_ub = np.vstack( A_list )
b_ub = np.hstack( b_list )

#regularity in time constraint
#A_list.append( partial_t )
#b_list.append( np.ones( partial_t.shape[0]) )
#A_list.append( -partial_t )
#b_list.append( np.ones( partial_t.shape[0]) )
#
##regularity in space constraint
#A_list.append( grad_k )
#b_list.append( 100*np.ones( grad_k.shape[0]))
#A_list.append( -grad_k )
#b_list.append( 100*np.ones( grad_k.shape[0]))

#boundary condition
#A_list = []
#b_list = []
#A_list.append( np.hstack( [np.zeros( (boundary_condition.shape[0] , Nx)), boundary_condition]))
#b_list.append( np.zeros( boundary_condition.shape[0]) )
#A_eq = np.vstack( A_list )
#b_eq = np.hstack( b_list )


#result = linprog(c,A_ub=A_ub,b_ub=b_ub, A_eq=A_eq,b_eq=b_eq)
result = linprog(c,A_ub=A_ub,b_ub=b_ub) 

w = result.x[0:Nx]
v = result.x[Nx:Nx+Nx*Nt]

plt.plot( result.x)
plt.show()
v.resize((Nt,Nx))
if(result.success):
    plt.imshow(v,\
            cmap='Greys',\
            interpolation='nearest')
    plt.show()
    plt.plot(x,w)
    plt.show()
print result.message
print result.status
