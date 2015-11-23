from time import time
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
from scipy.integrate import ode
# FOR PLOTTING ANYTHING UNCOMMENT THE FOLLOWING LINE
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#THERE IS AN INDEXING MISMATCH BETWEEN HOW PSI IS CALCULATED AND OUR OPERATOR

# CONSTANTS 
# modes is a list of integer tuples
p = 7
modes = [ (k1,k2) for k1 in range(- 2**p , 2**p+1 ) for k2 in range(- 2**p , 2**p + 1)]
res = 2**(p+1)+1

print 'There are %d modes' %(len(modes))
print 'So our operator has %d entries' %(len(modes)**2)

# These functions allow us to do basic arithmetic with the modes
add_tuple = lambda m,n: tuple( x+y for (x,y) in zip(m,n) )
sub_tuple = lambda m,n: tuple( x-y for (x,y) in zip(m,n) )
dot = lambda m,n: sum( x*y for (x,y) in zip(m,n) )

# f_hat is a dictionary of complex valued multi-dimensional fourier coefficients
# of the form {k:coefficient} where k is a frequency (possibly multi-dimensional)
# and coefficient is the corresponding Fourier coefficient (also multi-dimensional)
B = 1.0
C = 1.0
f_hat = {}
f_hat[(0,1)] = C*0.5
f_hat[(0,-1)] = C*0.5

g_hat = {}
g_hat[(1,0)] = B*0.5j
g_hat[(-1,0)] = -B*0.5j

# vector field is X(x) = f(x) d/dx
# f_hat is the Fourier transform of f

# Arrays for plotting stuff
temp = np.linspace(-0.5,0.5,res)
X_grid,Y_grid = np.meshgrid(temp,temp)

def ode_func(t,y):
    # This is for initializing our integration.
    global Operator
    return Operator.dot(y)
    #return np.zeros( y.shape)

def ode_jac(t,y):
    # we are not useing this function at the moment
    # although it should be supplied to ode() for optimal performance
    global Operator
    return Operator

def update(*args ):
    global Koopman_gen,FP_gen,Hilbert_gen, psi_hat,im
    t0 = time()
    psi_hat += 0.01*Hilbert_gen.dot(psi_hat)
    t1 = time()
    print 'it took %f seconds to produce psi_hat' %(t1-t0)
    psi = IFT(psi_hat)
    t2 = time()
    print 'it took %f seconds to produce psi' %(t2-t1)
    im.set_array( psi.real**2 )
    t2 = time()
    print 'it took %f seconds to plot' %(t2-t1)
    return im,

def conv_basis(k):
    # Takes in a mode and outputs a convolution matrix associated
    # to the activation of that mode.
    global res
    e_0 = sparse.eye(res,k=k[0])
    e_1 = sparse.eye(res,k=k[1])
    return sparse.kron(e_0,e_1)


def get_conv_mat( g_hat ):
    # outputs a convolution operator associated to 1D array g_hat
    out = sparse.dia_matrix((res**2,res**2),dtype=complex)
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
    sp_id = sparse.eye(res)
    if dimension_of_translation == 0:
        out = sparse.kron(ddx,sp_id)
    else:
        out = sparse.kron(sp_id,ddx)
    return out.todia()

def get_Koopman_generator():
    global modes,f_hat,g_hat
    N = len(modes)
    Op = sparse.dia_matrix( (N,N) , dtype = complex)
    ddx = get_translation_generator(0)
    f_cnv = get_conv_mat(f_hat)
    Op += f_cnv*ddx

    ddx = get_translation_generator(1)
    f_cnv = get_conv_mat(g_hat)
    Op += f_cnv*ddx
    return Op.todia()

def initialize_wave_function():
    # generates a wave function using fft
    global modes
    N_nodes = 2**(p+1) + 1 # should be equal to the number of Fourier modes
    x_arr = np.linspace(0.0,1.0,N_nodes)
    X,Y = np.meshgrid( x_arr , x_arr)
    sigma = 0.2
    psi = np.zeros(X.shape)
    # here we construct a wrapped Guassian
    x_bar = 0.5
    y_bar = 0.5
    for k1 in range(-7,7):
        for k2 in range(-7,7):
            if True or k1**2 + k2**2 < 7**2:
                s1 = k1*sigma
                s2 = k2*sigma
                psi += np.exp(- ((X-x_bar-k1)**2 + (Y-y_bar -k2)**2)/(2*sigma**2))
    plt.imshow(psi)
    plt.title('psi at t=0')
    plt.colorbar()
    plt.show()
    psi_hat = np.fft.fftn(psi)
    return psi_hat.flatten()

#Lets make some operators
print 'Let\'s make some operators'
t0=time()
Koopman_gen = get_Koopman_generator()
FP_gen = Koopman_gen.transpose().conj().copy()
Hilbert_gen = 0.5*Koopman_gen -0.5*FP_gen
print 'Done.  That took me %f seconds.' %(time()-t0)

nnz = Hilbert_gen.nnz
sparsity = (100.0 * nnz) / ( len(modes)**2)
print 'There are %d nonzero entries, sparsity = %f %%' %(nnz,sparsity)

psi_hat_initial = initialize_wave_function()
print 'integrating'
print 'p=%d' %p
Operator = Hilbert_gen
psi = np.fft.ifftn(psi_hat_initial.reshape(res,res))
dpsi_hat = Operator.dot(psi_hat_initial)
dpsi = np.fft.ifftn(dpsi_hat.reshape(res,res))
fig = plt.figure()
plt.imshow(dpsi.imag)
plt.colorbar()
plt.title('dpsi')
plt.show()
print np.max(np.abs(dpsi.imag) )
print np.max(np.abs(dpsi.real) )
print "The first number should be much smaller than the second"
ode_instance = ode( ode_func )
ode_instance.set_integrator('zvode',method='adams')
t0 = 0.0
ode_instance.set_initial_value(psi_hat_initial,t0)
#ode_instance.set_f_params(Operator)
#ode_instance.set_jac_params(Operator)
t1 = 1.0
dt = t1/90.
N_nodes = 2**(p+1)+1
fig_num = 0
while ode_instance.successful() and ode_instance.t < t1:
    ode_instance.integrate( ode_instance.t + dt )
    print 'Integrated up to t = %f' %(ode_instance.t)
    psi_hat = ode_instance.y
    psi = np.fft.ifftn( psi_hat.reshape(N_nodes,N_nodes))
    plt.imshow( np.abs(psi)**2,
            cmap = 'Greys',
            interpolation='bicubic',
            extent = [-0.5,0.5,-0.5,0.5])
    plt.grid(True)
    plt.title('t=%.2f' %ode_instance.t)
    plt.savefig('figures/figure_%d.png'%(fig_num))
    plt.clf()
    fig_num += 1
    print 'min( rho ) = %f' %(np.min(np.min( np.abs(psi)**2)))
    print 'total mass = %f' %(np.sum(np.sum( np.abs(psi)**2)))
#
#fig = plt.figure()
#im = plt.imshow( psi.real )
#ani = animation.FuncAnimation(fig, update, frames=30,interval=20,blit=True)
#ani.save('test.mp4')
