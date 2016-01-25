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
p = 2
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
f_hat,g_hat = {},{}

# THIS IS THE HAMILTONIAN COMPONENT
f_hat[( 0, 1)] = C/2
f_hat[( 0,-1)] = C/2
g_hat[( 1, 0)] = -B/2j
g_hat[(-1, 0)] = B/2j

# THIS IS THE DISSAPATIVE COMPONENT
D = 0.2
f_hat[( 1, 0)] = -D*B /2j
f_hat[(-1, 0)] = D*B/2j
g_hat[( 0, 1)] = -D*C/2
g_hat[( 0,-1)] = -D*C/2


def CB_flow( state , t ):
    global B,C,D
    x,y = state_to_XY( state)
    dx = C * np.cos(2*np.pi*y) - D*B*np.sin(2*np.pi*x)
    dy = -B * np.sin(2*np.pi*x) - D*C*np.cos(2*np.pi*y)
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

# Arrays for plotting stuff
temp = np.linspace(-0.5,0.5,res)
X_grid,Y_grid = np.meshgrid(temp,temp)
def initialize_wave_function():
    # generates a wave function using fft
    global modes
    N_nodes = 2**(p+1) + 1 # should be equal to the number of Fourier modes
    x_arr = np.linspace(0.0,1.0,N_nodes)
    X,Y = np.meshgrid( x_arr , x_arr)
    sigma_x = 0.3
    sigma_y = 0.3
    psi = np.zeros(X.shape)
    #here we construct a wrapped Guassian
    x_bar = 0.4
    y_bar = 0.4
    for k1 in range(-7,8):
        for k2 in range(-7,8):
            if True or k1**2 + k2**2 < 7**2:
                s1 = k1*sigma_x
                s2 = k2*sigma_y
                psi += np.exp( -(X-x_bar-k1)**2 / (2*sigma_x**2) - (Y-y_bar -k2)**2 /(2*sigma_y**2))
    psi_hat = np.fft.fftn(psi)
    N_samples = (100)**2
    rt2 = np.sqrt(2)
    x = ( np.random.randn(N_samples)*sigma_x/2 + x_bar ) 
    y = ( np.random.randn(N_samples)*sigma_y/2 + y_bar ) 
    psi = np.fft.ifftn(psi_hat)
    rho = psi**2
    rho_hat = np.fft.fftn( rho )
    plt.subplot(1,2,1)
    plt.imshow( np.abs(psi**2) , cmap='Greys' ,origin='lower', extent=[0,1,0,1] )
    plt.title('rho at t=0')
    plt.subplot(1,2,2,aspect='equal')
    plt.axis([0,1,0,1])
    plt.scatter( x%1,y%1 ,alpha=0.05, c= 'k', s=3. )
    plt.show()
    return psi_hat.flatten(), rho_hat.flatten(), x%1,y%1

#GET INITIAL CONDITIONS
psi_hat_initial,rho_hat_initial,x,y = initialize_wave_function()
print 'integrating'
print 'p=%d' %p
t0 = 0.0
t1 = 1.0
N_timesteps = 5
t = np.linspace(0.0 , t1 , N_timesteps)

#SOLVE USING MONTE-CARLO
state_0 = np.concatenate( (x , y) )
states = odeint( CB_flow ,state_0 , t )

#SOLVE USING HILBERT METHOD
Operator = -get_Hilbert_Op() 
ode_instance = ode( lambda t,y: Operator.dot(y) )
ode_instance.set_integrator('zvode',method='adams')
ode_instance.set_initial_value(psi_hat_initial,t0)
N_nodes = 2**(p+1)+1
psi_list = [np.fft.ifftn(psi_hat_initial.reshape(N_nodes,N_nodes)) ]
for i in range(1,N_timesteps):
    ode_instance.integrate( t[i] )
    print 'Integrated up to t = %f' %(ode_instance.t)
    psi_hat = ode_instance.y
    psi_list.append( np.fft.ifftn( psi_hat.reshape(N_nodes,N_nodes)) )

#SOLVE USING FP-CALLOCATION METHOD
Operator = -get_FP_Op() 
ode_instance = ode( lambda t,y: Operator.dot(y) )
ode_instance.set_integrator('zvode',method='adams')
ode_instance.set_initial_value(rho_hat_initial,t0)
N_nodes = 2**(p+1)+1
rho_list = [np.fft.ifftn(rho_hat_initial.reshape(N_nodes,N_nodes)) ]
for i in range(1,N_timesteps):
    ode_instance.integrate( t[i] )
    print 'Integrated up to t = %f' %(ode_instance.t)
    rho_hat = ode_instance.y
    rho_list.append( np.fft.ifftn( rho_hat.reshape(N_nodes,N_nodes)) )


#NOW TO MAKE A VIDEO
make_video = False
if make_video:
    f,ax = plt.subplots(1,3)
    for i in range(N_timesteps):
        #PLOT GN SPEC
        psi = psi_list[i]
        rho = rho_list[i]
        ax[0].imshow( np.abs(psi**2) ,
                cmap = 'Greys',
                origin='lower',
                interpolation='nearest',
                extent=[0,1,0,1])
        ax[0].grid(True)
        ax[0].set_title('GN spectral')    
        #PLOT MONTE CARLO
        x,y = state_to_XY( states[i] )
        ax[1].scatter( x%1,y%1 ,alpha=0.02, s=3.0, c = 'k')
        ax[1].axis([0,1,0,1])
        ax[1].set_aspect('equal')
        ax[1].grid(True)
        ax[1].set_title('Monte Carlo')

        #PLOT STANDARD SPEC
        ax[2].imshow( rho.real ,
                cmap = 'Greys',
                origin='lower',
                interpolation='nearest',
                extent=[0,1,0,1])
        ax[2].grid(True)
        ax[2].set_title('Standard spectral')
        #REMOVE TICKLABELS
        for k in range(3):
            ax[k].xaxis.set_ticklabels([]) 
            ax[k].yaxis.set_ticklabels([]) 
        plt.tight_layout()
        f.suptitle('t=%.2f'%t[i])
        plt.savefig('figures/figure_%d.png'%(i))
        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        print '%d of %d frames complete' %(i+1,N_timesteps)
        quit()

#NOW ON TO PLOTTING
if N_timesteps > 5:
    print "N timesteps is too large, quitting"
    quit()

f,ax = plt.subplots( 3, N_timesteps)
#f.set_size_inches(16,4)
#f.text( 0.02, 0.5, 't=%.2f'%t[i], fontsize=16)
for i in range(N_timesteps):
    #PLOT GN SPEC
    psi = psi_list[i]
    rho = rho_list[i]
    ax[0,i].imshow( np.abs(psi**2) ,
            cmap = 'Greys',
            origin='lower',
            interpolation='nearest',
            extent=[0,1,0,1])
    ax[0,i].grid(True)
    ax[0,i].set_title('t=%.2f'%t[i])    
    #PLOT MONTE CARLO
    x,y = state_to_XY( states[i] )
    ax[1,i].scatter( x%1,y%1 ,alpha=0.01, s=2.0, c = 'k')
    ax[1,i].axis([0,1,0,1])
    ax[1,i].set_aspect('equal')
    ax[1,i].grid(True)

    #PLOT STANDARD SPEC
    ax[2,i].imshow( rho.real ,
            cmap = 'Greys',
            origin='lower',
            interpolation='nearest',
            extent=[0,1,0,1])
    ax[2,i].grid(True)
    
    #REMOVE TICKLABELS
    for k in range(3):
        ax[k,i].xaxis.set_ticklabels([]) 
        ax[k,i].yaxis.set_ticklabels([]) 
    #plt.savefig('figures/figure_%d.png'%(fig_num))
    print 'min( rho ) = %f' %(np.min(np.min( np.abs(psi)**2)))
    print 'total mass = %f' %(np.sum(np.sum( np.abs(psi)**2)))

plt.tight_layout()
plt.show()
