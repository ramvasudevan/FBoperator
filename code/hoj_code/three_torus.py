from time import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
import scipy.io
from scipy.integrate import ode
from scipy.integrate import odeint
import random
# FOR PLOTTING ANYTHING UNCOMMENT THE FOLLOWING LINE
#from mpl_toolkits.mplot3d import axes3d
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
B = 0.5
C = 0.2

def ABC_flow( arr ,t):
    global A,B,C
    N = len(arr)
    x = arr[0:N/3]
    y = arr[N/3:2*N/3]
    z = arr[2*N/3:N]
    dx = A*np.sin(2*np.pi*z) + C*np.cos(2*np.pi*y)
    dy = B*np.sin(2*np.pi*x) + A*np.cos(2*np.pi*z)
    dz = C*np.sin(2*np.pi*y) + B*np.cos(2*np.pi*x)
    return np.hstack([dx,dy,dz])

#(u,v,w) = ( Asin(z)+Ccos(y), B sin(x) + A cos(z) , C sin(y) + B cos(x) )
f_hat = [{},{},{}]

f_hat[0][(0, 1,0)] = C*0.5
f_hat[0][(0,-1,0)] = C*0.5
f_hat[0][(0,0, 1)] = -A*0.5j
f_hat[0][(0,0,-1)] = A*0.5j
f_hat[1][( 1,0,0)] = -B*0.5j
f_hat[1][(-1,0,0)] = B*0.5j
f_hat[1][(0,0, 1)] = A*0.5
f_hat[1][(0,0,-1)] = A*0.5
f_hat[2][(0, 1,0)] = -C*0.5j
f_hat[2][(0,-1,0)] = C*0.5j
f_hat[2][( 1,0,0)] = B*0.5
f_hat[2][(-1,0,0)] = B*0.5

# Arrays for plotting stuff
temp = np.linspace(-0.5,0.5,res)
X_grid,Y_grid,Z_grid= np.meshgrid(temp,temp,temp)

def conv_basis(k):
    # Takes in a mode and outputs a convolution matrix associated
    # to the activation of that mode.
    global res
    data = np.ones([2,res])
    offsets = np.array( [2**p+1, -2**p])
    Shift = sparse.dia_matrix( (data,offsets) , shape=(res,res),dtype=complex )
    iShift = Shift.transpose()
    out = Shift.dot( sparse.eye(res,k=k[0]).dot( iShift) )
    for d in range(1,DIM):
        e =  Shift.dot( sparse.eye(res,k=k[d]).dot( iShift) )
        store = sparse.kron(out,e)
        out = store.copy()
    return out

def get_conv_mat( g_hat ):
    # outputs a the (conjugation of) the 
    # convolution operator associated to 1D array g_hat
    global modes
    N = len(modes)
    out = sparse.dia_matrix((N,N),dtype=complex)
    for k in g_hat.keys():
        e_k = conv_basis(k)
        out += g_hat[k].conjugate()*e_k
    return out

def get_translation_generator( dimension_of_translation ):
    # creates the infinitesimal generator along the ath coordinate
    global modes,p,res
    freq = np.fft.fftfreq(res)
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

ddx = []
f_cnv = []
for d in range(3):
    ddx.append( get_translation_generator(d) )
    f_cnv.append( get_conv_mat(f_hat[d]) )

def get_Koopman_Op():
    def mv(v):
        global f_cnv,ddx
        w = np.zeros(v.size,dtype=complex)
        for d in range(3):
            w += f_cnv[d].dot( ddx[d].dot(v) )
        return w
    return LinearOperator( f_cnv[0].shape , matvec=mv )

def get_Frobenius_Perron_Op():
    def mv(v):
        global f_cnv,ddx
        w = np.zeros(v.size,dtype=complex)
        for d in range(3):
            w += ddx[d].dot( f_cnv[d].dot(v) )
        return w
    return LinearOperator( f_cnv[0].shape , matvec=mv )

def get_Hilbert_Op():
    def mv(v):
        global f_cnv,ddx
        w1 = np.zeros(v.size,dtype=complex)
        w2 = np.zeros(v.size,dtype=complex)
        for d in range(3):
            w1 += f_cnv[d].dot( ddx[d].dot(v) )
            w2 += ddx[d].dot( f_cnv[d].dot(v) )
        return 0.5*w1+0.5*w2
    return LinearOperator( f_cnv[0].shape , matvec=mv )

def Gaussian_in_Fourier(particles=False):
    #generates the Fourier transform of a Guassian
    global modes,res
    N = len(modes)
    sigma = [ 0.2 , 0.3 , 0.3 ]
    modes_arr = np.array(modes)
    freq = np.fft.fftfreq(res)
    unos = np.ones(res)
    fx = np.einsum('i,j,k->ijk',freq,unos,unos)
    fy = np.einsum('i,j,k->ijk',unos,freq,unos)
    fz = np.einsum('i,j,k->ijk',unos,unos,freq)
    r2 = (fx*sigma[0])**2 + (fy*sigma[1])**2 + (fz*sigma[2])**2
    out = np.exp( -2*np.pi**2 * r2 * (res**2) )
    print modes_arr.shape
    if particles:
        N_pts = 15**3
        X = sigma[0]*np.random.randn(N_pts) / 2#np.sqrt(2)
        Y = sigma[1]*np.random.randn(N_pts) / 2#np.sqrt(2)
        Z = sigma[2]*np.random.randn(N_pts) / 2#np.sqrt(2)
        return out,X%1,Y%1,Z%1
    return out

#Lets make some operators
print 'Let\'s make some operators'
t0=time()
Operator = get_Hilbert_Op()
#Operator = ddx[0]
print 'Done.  That took me %f seconds.' %(time()-t0)
print 'initializing psi with p=%d and res=%d' %(p,res)
t0 = time()
psi_hat,X_pts,Y_pts,Z_pts = Gaussian_in_Fourier(particles=True)
print 'Done.  That took me %f seconds.' %(time()-t0)
print 'Integrating'
psi = np.fft.ifftn(psi_hat.reshape(res,res,res))

#LET'S DISPLAY THE INITIAL CONDITION BEFORE PROCEEDING
f,ax = plt.subplots(1,2)
ax[0].imshow( (np.abs(psi)**2).sum(2).transpose() ,
        cmap='Greys',
        extent=[0,1,0,1],
        interpolation='nearest',
        origin = 'lower')
ax[0].grid(True)
ax[0].set_title('GN spectral')
ax[1].scatter(X_pts%1,Y_pts%1,alpha=0.15,s=3,c='k')
ax[1].set_aspect('equal')
ax[1].axis([0,1,0,1])
ax[1].set_title('Monte Carlo')
ax[1].grid(True)
plt.show()
f.clear()

N_timesteps = 30
t = np.linspace(0,1,N_timesteps+1 )
r = ode( lambda t,y: -Operator.dot(y) ).set_integrator('zvode', method='bdf')
r.set_initial_value( psi_hat.flatten() , t[0] )


#INTEGRATING POINT CLOUD
states = odeint( ABC_flow , np.hstack( [X_pts,Y_pts,Z_pts] ) , t )
N_pts = states[0].size / 3
f,ax = plt.subplots(1,2)
for i in range(N_timesteps):
    psi = np.fft.ifftn( r.y.reshape(res,res,res) )
    r.integrate(t[i+1])
    if r.successful():
        print 'success'
    else:
        print 'fail'
    ax[0].imshow( (np.abs(psi)**2).sum(2).transpose() ,
            cmap = 'Greys',
            interpolation='nearest',
            extent = [ 0. , 1. , 0. , 1. ],
            origin = 'lower')
    ax[0].grid(True)
    ax[0].set_title('GN Spectral')
    X_pts = states[i,0:N_pts]
    Y_pts = states[i,N_pts:2*N_pts]
    ax[1].scatter( X_pts % 1 , Y_pts % 1, alpha = .12, s=3., c='k')
    ax[1].axis([0,1,0,1])
    ax[1].set_aspect('equal')
    ax[1].set_title('Monte Carlo')
    ax[1].grid(True)
    plt.savefig('./figures/figure_{:,d}'.format(i))
    ax[0].clear()
    ax[1].clear()
    print "print frame %d of %d" %(i,N_timesteps)
