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
p = 5
modes = [ (k1,k2) for k1 in range(- 2**p , 2**p+1 ) for k2 in range(- 2**p , 2**p + 1)]

# These functions allow us to do basic arithmetic with the modes
add_tuple = lambda m,n: tuple( x+y for (x,y) in zip(m,n) )
sub_tuple = lambda m,n: tuple( x-y for (x,y) in zip(m,n) )
dot = lambda m,n: sum( x*y for (x,y) in zip(m,n) )

# f_hat is a dictionary of complex valued multi-dimensional fourier coefficients
# of the form {k:coefficient} where k is a frequency (possibly multi-dimensional)
# and coefficient is the corresponding Fourier coefficient (also multi-dimensional)
f_hat = {}
f_hat[(0,0)] = np.array( [0.0 , 0.0 ] , dtype=complex)
f_hat[(1,0)] = np.array( [ -0.5j , 0] )
f_hat[(-1,0)] = np.array( [ 0.5j , 0] )
f_hat[(0,1)] = np.array( [ 0.0 , -0.5j] )
f_hat[(0,-1)] = np.array( [ 0.0 , 0.5j] )

# vector field is X(x) = f(x) d/dx
# f_hat is the Fourier transform of f

# Arrays for plotting stuff
temp = np.linspace(-0.5,0.5,100)
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

def get_f_conv(a):
    global f_hat
    N = len(modes)
    offsets = ()
    data = np.zeros([len(f_hat) , N],dtype=complex)
    for k in f_hat.keys():
        offsets += (modes.index(k),)
        data[k,:] = f_hat[k][a]*np.ones(N)
        
    plt.matshow( data.imag )
    plt.show()
    out= sparse.dia_matrix( (data, offsets) \
            , shape=(N,N),dtype=complex)
    plt.matshow( out.todense().imag )
    plt.show()
    return out

def get_translation_generator( dimension_of_translation ):
    # creates the infinitesimal generator along the ath coordinate
    global modes,p
    #data = np.zeros(len(modes),dtype=complex)
    #store = np.zeros( [2**(p+1) + 1 , 2**(p+1)+1],dtype=complex )
    #offsets = np.array([0])
    #for k1 in range(-2**p,2**p+1):
    #    for k2 in range(-2**p,2**p+1):
    #        store[k1,k2] = 2*np.pi*1j*k1
    #data = store.flatten()
    #for k in modes:
    #    store[k[0],k[1]] = 2*np.pi*1j*k[a]
    #data = store.flatten()
    res = 2**(p+1) + 1
    freq = np.fft.fftfreq(res)
    ddx = sparse.dia_matrix(( 2*np.pi*1j*freq ,
        np.array([0])),
        shape=(res,res) , dtype=complex)
    sp_id = sparse.eye(res)
    if dimension_of_translation == 0:
        out = sparse.kron(ddx,sp_id)
    else:
        out = sparse.kron(sp_id,ddx)
    #return sparse.dia_matrix((data,offsets), shape=(len(modes),len(modes)),dtype=complex)
    return out.todia()

def get_Koopman_generator():
    # SOMETHING IS REALLY WRONG HERE I THINK
    N = len(modes)
    Op = sparse.dia_matrix( (N,N) , dtype = complex)
    for a in range(0,2):
        ddx = get_translation_generator(a)
        f_cnv = get_f_conv(a)
        Op += 0.5*f_cnv*ddx
    return Op

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
#print 'Let\'s make some operators'
#t0=time()
#Koopman_gen = get_Koopman_generator()
#FP_gen = ((Koopman_gen.transpose()).copy()).conj()
#Hilbert_gen = 0.5*Koopman_gen -0.5*FP_gen
#print 'Done.  That took me %f seconds.' %(time()-t0)

psi_hat_initial = initialize_wave_function()
print 'integrating'
t = time()
t_arr = np.linspace(0,1,30)
dimension = 1
Operator= get_translation_generator(dimension)
#Operator = Koopman_gen

print 'Let us test this translation operator'
i = (1,-3)
j = (1,-3)
print 'i = (%d,%d)' %(i)
print '2*pi*i[%d]*1j = %f * 1j' %(dimension,2*np.pi*i[dimension])
print 'j = (%d,%d)' %(j)
print '2*pi*j[%d]*1j = %f * 1j' %(dimension,2*np.pi*j[dimension])
print Operator.todense()[ modes.index(i) , modes.index(j) ]

print 'Is our operator anti-Hermetian?'
A = Operator.todense()
print 'It is if %f is small.'\
        %(np.max(np.max( np.abs(A.conj() + A.transpose() ))))

print 'p=%d' %p
res = 2**(p+1)+1
psi = np.fft.ifftn(psi_hat_initial.reshape(res,res))
dpsi_hat = Operator.dot(psi_hat_initial)
dpsi = np.fft.ifftn(dpsi_hat.reshape(res,res))
fig = plt.figure()
plt.imshow(dpsi.real)
plt.colorbar()
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
t1 = 9.0
dt = 1.5
N_nodes = 2**(p+1)+1
while ode_instance.successful() and ode_instance.t < t1:
    ode_instance.integrate( ode_instance.t + dt )
    print 'Integrated up to t = %f' %(ode_instance.t)
    psi_hat = ode_instance.y
    psi = np.fft.ifftn( psi_hat.reshape(N_nodes,N_nodes))
    plt.imshow(psi.real)
    plt.show()
print 'Done.  That took %f seconds.' %(time()-t)
#
#fig = plt.figure()
#im = plt.imshow( psi.real )
#ani = animation.FuncAnimation(fig, update, frames=30,interval=20,blit=True)
#ani.save('test.mp4')
