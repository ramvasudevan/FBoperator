from time import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import linalg
from scipy.integrate import ode
from scipy.integrate import odeint
# FOR PLOTTING ANYTHING UNCOMMENT THE FOLLOWING LINE
import matplotlib.pyplot as plt


# CONSTANTS 
# modes is a list of integer tuples
p = 6
modes = [ (k1,k2) for k1 in range(- 2**p , 2**p+1 ) for k2 in range(- 2**p , 2**p + 1)]
res = 2**(p+1)+1

print 'There are %d modes' %(len(modes))
print 'So our operator has %d entries' %(len(modes)**2)

#CONVENTIONS:  DFT[ a ]_k = \sum_{x} a(x) exp( - 2*pi*i*x*k)
#So we must divide by n in order for this sum to approximate an integral.
# IDFT[ A ]_x = \frac{1}{n} \sum_{k} A_k exp( 2*pi*i*x*k)
#So we must multiply by n to be consistent with the above.

# f_hat and g_hat are dictionaries of fourier coefficients
# of the form {k:coefficient} where k is a (multi)frequency
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

def CB_flow( state , t ):
    global B,C
    x,y = state_to_XY( state)
    dx = C * np.cos(2*np.pi*y)
    dy = -B * np.sin(2*np.pi*x)
    return np.concatenate((dx,dy))

def translate_right( state , t):
    x,y = state_to_XY(state)
    return np.concatenate( (np.ones(x.size) , np.zeros(y.size)) )
    #return np.concatenate( (np.zeros(x.size) , np.ones(y.size)) )

def state_to_XY( state):
    x = state[0:len(state)/2]
    y = state[len(state)/2:]
    return x,y

def conv_basis(k):
    # Takes in a mode and outputs a convolution matrix associated
    # to the activation of that mode.
    global p,res
    data = np.ones([2,res])
    offsets = np.array( [2**p+1, -2**p])
    Shift = sparse.dia_matrix( (data,offsets) , shape=(res,res) )
    iShift = Shift.transpose()
    e_0 = Shift.dot( sparse.eye(res,k=k[0]).dot( iShift ) )
    e_1 = Shift.dot( sparse.eye(res,k=k[1]).dot( iShift ) )
    return sparse.kron(e_1,e_0)

def get_conv_mat( g_hat ):
    # outputs a convolution operator associated to 1D array g_hat
    out = sparse.dia_matrix((res**2,res**2),dtype=complex)
    for k in g_hat.keys():
        e_k = conv_basis(k)
        out += g_hat[k].conjugate()*e_k
    return out

def ddx_1d( dimension_of_translation ):
    # creates the infinitesimal generator along the ath coordinate
    global modes,p,res
    freq = np.fft.fftfreq(res)
    #I'm not sure why the following is multiplied by res
    ddx =  res*sparse.dia_matrix(( 2*np.pi*1j*freq ,
        np.array([0])),
        shape=(res,res) , dtype=complex)
    sp_id = sparse.eye(res)
    if dimension_of_translation == 1: 
        out = sparse.kron(ddx,sp_id)
    else:
        out = sparse.kron(sp_id,ddx)
    return out.todia()

f_cnv = get_conv_mat(f_hat)
g_cnv = get_conv_mat(g_hat)
ddx = ddx_1d(0)
ddy = ddx_1d(1)

def get_Koopman_Op():
    global modes,f_cnv,g_cnv,ddx,ddy
    N = len(modes)
    def mv(v):
        w1 =  f_cnv.dot( ddx.dot(v) )
        w2 = g_cnv.dot(ddy.dot(v))
        return w1+w2
    return LinearOperator( (N,N) , matvec=mv)

def get_FP_Op():
    global modes,f_cnv,g_cnv,ddx,ddy
    N = len(modes)
    def mv(v):
        w1 =  ddx.dot( f_cnv.dot(v) )
        w2 = ddy.dot( g_cnv.dot(v))
        return w1+w2
    return LinearOperator( (N,N) , matvec=mv)

def get_Hilbert_Op():
    global modes,f_cnv,g_cnv,ddx,ddy
    N=len(modes)
    def mv(v):
        w1 =  f_cnv.dot( ddx.dot(v) )
        w2 = g_cnv.dot(ddy.dot(v))
        w3 =  ddx.dot( f_cnv.dot(v) )
        w4 = ddy.dot( g_cnv.dot(v))
        return 0.5*(w1+w2+w3+w4)
    return LinearOperator( (N,N) , matvec=mv)

def initialize_wave_function():
    # generates a wave function using fft
    global modes
    N_nodes = 2**(p+1) + 1 # should be equal to the number of Fourier modes
    x_arr = np.linspace(0.0,1.0,N_nodes)
    X,Y = np.meshgrid( x_arr , x_arr)
    sigma_x = 0.3
    sigma_y = 0.1
    psi = np.zeros(X.shape)
    # here we construct a wrapped Guassian
    x_bar = 0.2
    y_bar = 0.9
    for k1 in range(-7,8):
        for k2 in range(-7,8):
            if True or k1**2 + k2**2 < 7**2:
                s1 = k1*sigma_x
                s2 = k2*sigma_y
                psi += np.exp( -(X-x_bar-k1)**2 / (2*sigma_x**2) - (Y-y_bar -k2)**2 /(2*sigma_y**2))
    psi_hat = np.fft.fftn(psi)
    N_samples = (50)**2
    x = ( np.random.randn(N_samples)*sigma_x/2 + x_bar ) 
    y = ( np.random.randn(N_samples)*sigma_y/2 + y_bar ) 
    psi = np.fft.ifftn(psi_hat)
    rho = (psi**2).real
    plt.subplot(1,2,1)
    plt.imshow( rho , cmap='Greys' ,origin='lower', extent=[0,1,0,1] )
    plt.title('rho at t=0')
    plt.subplot(1,2,2,aspect='equal')
    plt.axis([0,1,0,1])
    plt.scatter( x%1,y%1 ,alpha=0.035, c = 'k')
    plt.show()
    return psi_hat.flatten(),x%1,y%1

#Lets make some operators
psi_hat_initial,x,y = initialize_wave_function()
print 'integrating'
print 'p=%d' %p
#Operator = -ddx
Operator = -get_Hilbert_Op() 
ode_instance = ode( lambda t,y: Operator.dot(y) )
ode_instance.set_integrator('zvode',method='adams')
t0 = 0.0
t1 = 1.0
ode_instance.set_initial_value(psi_hat_initial,t0)
t_span = np.linspace(0.0 , t1 , 10)
dt = t_span[1]
state_0 = np.concatenate( (x , y) )
states = odeint( CB_flow ,state_0 , t_span )
#states = odeint( translate_right ,state_0 , t_span )
N_nodes = 2**(p+1)+1
fig_num = 0
x_span = np.linspace( 0 , 1. , N_nodes)
X_mesh,Y_mesh = np.meshgrid( x_span,x_span)
while ode_instance.successful() and ode_instance.t < t1:
    ode_instance.integrate( ode_instance.t + dt )
    print 'Integrated up to t = %f' %(ode_instance.t)
    psi_hat = ode_instance.y
    psi = np.fft.ifftn( psi_hat.reshape(N_nodes,N_nodes))
    plt.subplot(1,2,1)
    plt.imshow( np.abs(psi)**2,
            cmap = 'Greys',
            origin='lower',
            interpolation='nearest',
            extent=[0,1,0,1])
    plt.grid(True)
    plt.title('t=%.2f' %ode_instance.t)
    fig_num += 1
    x,y = state_to_XY( states[fig_num] )
    plt.subplot(1,2,2,aspect='equal')
    plt.axis([0,1,0,1])
    plt.scatter( x%1,y%1 ,alpha=0.035, c = 'k')
    plt.show()
    #plt.savefig('figures/figure_%d.png'%(fig_num))
    plt.clf()
    print 'min( rho ) = %f' %(np.min(np.min( np.abs(psi)**2)))
    print 'total mass = %f' %(np.sum(np.sum( np.abs(psi)**2)))
