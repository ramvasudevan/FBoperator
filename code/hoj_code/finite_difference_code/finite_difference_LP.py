import numpy as np
from scipy.optimize import linprog
from matplotlib import pyplot as plt

Nx = 12
Nt = 20
MIN_X = -1.0
MAX_X = 1.0
MIN_T = 0.0
MAX_T = 1.0

# When choosing these parameters note the CFL condition.
# We want u*dt / dx <= 1, where u is the max speed of the vector-field.
# This means we want u*(Nx-1)/(Nt-1) <=1

span_x = (MIN_X,MAX_X,Nx)
span_t = (MIN_T,MAX_T,Nt)
x = np.linspace( *span_x )
t = np.linspace( *span_t )
dx = x[1]-x[0]
dt = t[1]-t[0]

def derivative_op(x_min,x_max,N):
    #produces a one-dimesional derivative operator
    #sends N-vectors to (N-1)-vectors
    dx = (x_max - x_min)/float(N-1)
    out = (np.eye(N-1,N,k=1)-np.eye(N-1,N))/dx
    return out

def midpoint_op(N):
    #produces a one-dimensional mid-point interpolant.
    #sends N-vectors to (N-1)-vectors
    return 0.5*(np.eye(N-1,N,k=1) + np.eye(N-1,N))

def higher_order_derivative_op(k, x_min, x_max, N):
    assert(k>0)
    out = np.eye(N)
    for j in range(k):
        out = derivative_op(x_min,x_max,N).dot(out)
        dx = (x_min - x_max) / float(N-1)
        x_min += dx/2
        x_max -= dx/2
        N = N-1
    return out

def unit_test_derivative_op(min_x,max_x,N):
    x = np.linspace(min_x,max_x,N)
    y = np.sin(x)
    ddx = derivative_op(min_x,max_x,N)
    dy_dx = np.dot(ddx,y)
    x_mid = midpoint_op(N).dot( x )
    plt.plot( x_mid , dy_dx-np.cos(x_mid) , 'r-')
    plt.title('error')
    plt.show()
    return 0

#unit_test_derivative_op( *span_x)
#quit()

from scipy.linalg import kron

partial_x = kron( midpoint_op(Nt), derivative_op( *span_x ) )
partial_t = kron( derivative_op( *span_t ), midpoint_op(Nx) )

def unit_test_partial_t():
    x = np.linspace( *span_x )
    x_mid = midpoint_op(Nx).dot(x)
    t = np.linspace( *span_t )
    t_mid = midpoint_op(Nt).dot(t)
    #consider the function x*(t**2)
    f = np.kron( t**2 , x )
    df_dt_computed = np.dot( partial_t , f )
    df_dt_computed.resize( Nt-1,Nx-1)
    df_dt = np.kron( 2*t_mid , x_mid )
    df_dt.resize(Nt-1,Nx-1)
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
    x_mid = midpoint_op(Nx).dot(x)
    t = np.linspace( *span_t )
    t_mid = midpoint_op(Nt).dot(t)
    #consider the function exp(t)*sin(x)
    f = np.kron( np.exp(t) , np.sin(2*np.pi*x) )
    df_dx_computed = np.dot( partial_x , f )
    df_dx_computed.resize( Nt-1,Nx-1)
    df_dx = np.kron( np.exp(t_mid) , 2*np.pi*np.cos(2*np.pi*x_mid) )
    df_dx.resize(Nt-1,Nx-1)
    plt.subplot(2,1,1)
    plt.imshow( df_dx )
    plt.ylabel('x')
    plt.subplot(2,1,2)
    plt.imshow( df_dx_computed )
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()
    return 0

L = partial_t + 0.5*partial_x

#Create restriction operators
#X_T = [-.5 , 0.5]
XT_indices = [ i for i in xrange(Nx) if (i > Nx/4) and (i < 3*Nx/4) ]
XTc_indices = [i for i in xrange(Nx) if i not in XT_indices ]

Ix = np.eye(Nx)
at_time_T = np.zeros( Nt )
at_time_T[Nt-1] = 1.0
restrict_XT = np.kron(at_time_T, Ix[ XT_indices,:])
restrict_XTc = np.kron(at_time_T,  Ix[XTc_indices,:] )

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
delta_T = np.zeros(Nt)
delta_T[Nt-1] = 1.0
int_X = np.ones(Nx)*dx
int_X[0] = dx/2
int_X[-1] = dx/2
Ix = np.eye(Nx)
#c(f) = int_X f(0,x) dx
c = np.kron( delta_0 , int_X )

A_list = []
b_list = []

#Louiville constraint, L[v] <= 0
A_list.append( L  )
b_list.append( np.zeros(L.shape[0]) )

#f(T,x) >= 1 for x in X_T
A_list.append( - np.kron( delta_T , Ix[XT_indices,:] ) )
b_list.append( - np.ones( len(XT_indices) ) )

#f(T,x) >= 0 on X_T complement at time T
A_list.append( - np.kron( delta_T , Ix[XTc_indices,:] ) )
b_list.append( - np.zeros( len(XTc_indices) ) )

#grad_k f <= 1
#A_list.append( grad_k )
#b_list.append( 100*np.ones(Nx*Nt) )
#A_list.append( -grad_k )
#b_list.append( 100*np.ones(Nx*Nt) )

A_ub = np.vstack( A_list )
b_ub = np.hstack( b_list )

A_list=[]
b_list=[]
#boundary condition
on_left = np.zeros(Nx)
on_left[0] = 1.0
boundary_condition = np.kron( np.eye(Nt-1,Nt,k=1) , on_left )
A_list.append( boundary_condition)
b_list.append( np.zeros(Nt-1) )
A_eq = np.vstack( A_list )
b_eq = np.hstack( b_list )

print A_eq.shape

result = linprog(c,A_ub=A_ub,b_ub=b_ub, A_eq=A_eq,b_eq=b_eq )
#result = linprog(c,A_ub=A_ub,b_ub=b_ub) 

print result.message
print "status = {:d}.".format(result.status)

if(result.success):
    v = result.x
    v.resize((Nt,Nx))
    plt.imshow(v,\
            cmap='Greys',\
            interpolation='nearest')
    plt.show()
