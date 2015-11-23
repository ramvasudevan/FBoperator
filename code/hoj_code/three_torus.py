from time import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import scipy.io
#from scipy.integrate import ode
import random
# FOR PLOTTING ANYTHING UNCOMMENT THE FOLLOWING LINE
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

# CONSTANTS 
# modes is a list of integer tuples
p = 5

res = 2**(p+1)+1
freq = np.fft.fftfreq(res)
modes = [ (f1,f2,f3)    for f1 in freq
                        for f2 in freq
                        for f3 in freq]

print 'There are %e modes' %(len(modes))
print 'So our operator has %e entries' %(len(modes)**2)

# These functions allow us to do basic arithmetic with the modes
add_tuple = lambda m,n: tuple( x+y for (x,y) in zip(m,n) )
sub_tuple = lambda m,n: tuple( x-y for (x,y) in zip(m,n) )
dot = lambda m,n: sum( x*y for (x,y) in zip(m,n) )

# f_hat is a dictionary of complex valued multi-dimensional fourier coefficients
# A > B > C > 0 with C << 1 yields chaos
DIM = 3
A = 1.0
B = 0.7
C = 0.5

def ABC_flow( arr ,t):
    global A,B,C
    x = arr[0]
    y = arr[1]
    z = arr[2]
    dx = A*np.sin(z) + C*np.cos(y)
    dy = B*np.sin(x) + A*np.cos(z)
    dz = C*np.sin(y) + B*np.cos(x)
    return np.array([dx,dy,dz])

#(u,v,w) = ( Asin(z)+Ccos(y), B sin(x) + A cos(z) , C sin(y) + B cos(x) )
f_hat = [{},{},{}]

#f_hat[0][(0,0,0)] = 1.0
f_hat[0][(0, 1,0)] = C*0.5
f_hat[0][(0,-1,0)] = C*0.5
f_hat[0][(0,0, 1)] = A*0.5j
f_hat[0][(0,0,-1)] = -A*0.5j
#f_hat[1][(0,0,0)] = 1.0
f_hat[1][( 1,0,0)] = B*0.5j
f_hat[1][(-1,0,0)] = -B*0.5j
f_hat[1][(0,0, 1)] = A*0.5
f_hat[1][(0,0,-1)] = A*0.5
f_hat[2][(0, 1,0)] = C*0.5j
f_hat[2][(0,-1,0)] = -C*0.5j
f_hat[2][( 1,0,0)] = B*0.5
f_hat[2][(-1,0,0)] = B*0.5

# Arrays for plotting stuff
temp = np.linspace(-0.5,0.5,res)
X_grid,Y_grid,Z_grid= np.meshgrid(temp,temp,temp)

def cayley( A , y ):
    #outputs cay(A).dot(y)
    global modes
    Id = sparse.eye( len(modes), dtype=complex)
    store = spsolve(Id-A,y)
    return A.dot(y)+y

def ode_func(t,y):
    # This is for initializing our integration.
    global Operator
    return Operator.dot(y)
    #return np.zeros( y.shape)

def integrate_point_cloud(n_points,t_span):
    # initializes a normally distributed point cloud and advects it by the ABC flow
    points = np.random.randn(n_points,3)*0.25 % 1.0
    out = np.zeros([ n_points , len(t_span), 3])
    for point_index in range(0,n_points):
        y_initial = points[point_index,:]
        #out[point_index] = ode
    return 0

def conv_basis(k):
    # Takes in a mode and outputs a convolution matrix associated
    # to the activation of that mode.
    global res
    out = sparse.eye(res,k=k[0])
    for d in range(1,DIM):
        e =  sparse.eye(res,k=k[d])
        store = sparse.kron(out,e)
        out = store.copy()
    return out

def get_conv_mat( g_hat ):
    # outputs a convolution operator associated to 1D array g_hat
    global modes
    N = len(modes)
    out = sparse.dia_matrix((N,N),dtype=complex)
    for k in g_hat.keys():
        e_k = conv_basis(k)
        out += g_hat[k]*e_k
    return out

def get_translation_generator( dimension_of_translation ):
    # creates the infinitesimal generator along the ath coordinate
    global modes,p,res
    freq = np.fft.fftfreq(res)
    #I'm not sure why the following is multiplied by res
    ddx = res*sparse.dia_matrix(( 2*np.pi*1j*freq ,
        np.array([0])),
        shape=(res,res) , dtype=complex)
    sp_id = sparse.eye(res, dtype=complex)
    if dimension_of_translation == 0:
        out = ddx
    else:
        out = sp_id
    for d in range(1,DIM):
        if dimension_of_translation == d:
            store = sparse.kron(out,ddx)
        else:
            store = sparse.kron(out,sp_id)
        out = store.copy()
    return out.todia()

def get_Koopman_generator():
    global modes,f_hat
    N = len(modes)
    Op = sparse.dia_matrix( (N,N) , dtype = complex)
    for d in range(0,DIM):
        ddx = get_translation_generator(d)
        f_cnv = get_conv_mat(f_hat[d])
        Op += f_cnv*ddx
    return Op.todia()

def initialize_wave_function():
    # generates a wave function using fft
    global modes, res
    N = len(modes)
    x_arr = np.linspace(0.0,1.0,res)
    X = np.meshgrid( x_arr , x_arr , x_arr)
    sigma = 0.75
    psi = np.zeros([res,res,res],dtype=complex)
    # here we construct a wrapped Guassian
    x_bar = np.array([0.5,0.5,0.5])
    k_max = np.ceil( 7*sigma )
    k_list = ( np.array([k1,k2,k3])
            for k1 in np.arange(-k_max,k_max)
            for k2 in np.arange(-k_max,k_max)
            for k3 in np.arange(-k_max,k_max) )
    for k in k_list:
        r2 = np.zeros(X[0].shape)
        for d in range(0,3):
            r2 += (X[d] - x_bar[d] - k[d])**2
        psi += np.exp( -r2 / (2.0*sigma**2) )
    print '%f should be approximately 1.0' % psi.max()
    psi_hat = np.fft.fftn(psi)
    return psi_hat.flatten()

def Gaussian_in_Fourier(sigma):
    #generates the Fourier transform of a Guassian
    global modes,res
    N = len(modes)
    modes_arr = np.array(modes)
    print modes_arr.shape
    return np.exp( -sigma*(res**2*modes_arr**2).sum(1))

#Lets make some operators
print 'Let\'s make some operators'
t0=time()
Koopman_gen = get_Koopman_generator()
FP_gen = Koopman_gen.transpose().conj().copy()
Hilbert_gen = 0.5*Koopman_gen -0.5*FP_gen
print 'Done.  That took me %f seconds.' %(time()-t0)

nnz = Hilbert_gen.nnz
sparsity = (100.0 * nnz) / ( len(modes)**2)
print 'There are %e nonzero entries, sparsity = %f %%' %(nnz,sparsity)

print 'initializing psi with p=%d and res=%d' %(p,res)
t0 = time()
#psi_hat_initial = initialize_wave_function()
psi_hat_initial = Gaussian_in_Fourier(0.25)
print 'Done.  That took me %f seconds.' %(time()-t0)
print 'Integrating'
Operator = Hilbert_gen
#Operator = get_translation_generator(1) 

psi = np.fft.ifftn(psi_hat_initial.reshape(res,res,res))
dpsi_hat = Operator.dot(psi_hat_initial)
dpsi = np.fft.ifftn(dpsi_hat.reshape(res,res,res))

plt.imshow( psi.real.sum(2) ,
            cmap = 'Greys',
            interpolation='nearest',
            extent = [ 0. , 1. , 0. , 1. ])
plt.grid(True)
plt.title("psi")
plt.show()
plt.imshow( dpsi.real.sum(2) ,
            cmap = 'Greys',
            interpolation='nearest',
            extent = [ 0. , 1. , 0. , 1. ])
plt.grid(True)
plt.title("dpsi")
plt.show()

y_current = psi_hat_initial.copy()
dt = 0.005
t = 0.0
T = 1.0
frame_count = 0
x_arr = np.linspace(0,1,res)
X,Y = np.meshgrid( x_arr , x_arr )
while t < T:
    psi = np.fft.ifftn( y_current.reshape(res,res,res) )
    k1 = Operator.dot(y_current)
    k2 = Operator.dot(y_current + dt*k1/2.)
    k3 = Operator.dot(y_current + dt*k2/2.)
    k4 = Operator.dot(y_current + dt*k3)
    y_next = y_current + dt*(k1+2*k2+2*k3+k4)/6.
    if t > frame_count*T / 10.0:
        rho = np.abs(psi)**2
        fname = './density_data/frame_{:,d}'.format(frame_count)
        np.save(fname,rho)
        frame_count += 1
        print 'saved at time {:,.2f}'.format(t)
    t += dt
    y_current = y_next.copy()
   

#ode_instance = ode( ode_func )
#ode_instance.set_integrator('zvode')
#t0 = 0.0
#ode_instance.set_initial_value(psi_hat_initial,t0)
#ode_instance.set_f_params(Operator)
#ode_instance.set_jac_params(Operator)
#t1 = 1.00
#dt = 0.01
#fig_num = 0
#fig = plt.figure()
#while ode_instance.successful() and ode_instance.t < t1:
#    ode_instance.integrate( ode_instance.t + dt )
#    print 'Integrated up to t = %f' %(ode_instance.t)
#    psi_hat = ode_instance.y
#    psi = np.fft.ifftn( psi_hat.reshape(res,res,res))
#    plt.imshow( (np.abs(psi)**2).sum(0) ,
#            cmap = 'Greys',
#            interpolation='nearest',
#            extent = [-0.5,0.5,-0.5,0.5])
#    plt.grid(True)
#    plt.title('projection onto yz-plane, t=%.2f' %(ode_instance.t))
#    plt.savefig('figures/figure_%d.png'%(fig_num))
#    fig_num += 1
#    print 'min( rho ) = %f' %(np.min(np.min( np.abs(psi)**2)))
#    print 'total mass = %f' %(np.sum(np.sum( np.abs(psi)**2)) / (res)**3 )
