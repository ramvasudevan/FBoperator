import numpy as np
from scipy import sparse
from scipy.integrate import ode
import matplotlib.pyplot as plt

N = 32

def get_ddx(N):
    diag = 1j*np.arange(-N,N+1)
    data = np.zeros([1,2*N+1],dtype=complex)
    data[0,:] = diag
    offsets = np.zeros(1)
    return sparse.dia_matrix((data,offsets),shape=(2*N+1,2*N+1),dtype=complex)

def mult_by_e( k ):
    return sparse.eye(2*N+1,k=-k,dtype=complex, format='dia')

def mult_by_sin(k):
    return (mult_by_e(k) - mult_by_e(-k))/(2j)

def mult_by_cos(k):
    return (mult_by_e(k) + mult_by_e(-k))/2.

x = np.linspace(-np.pi, np.pi,200)
def ift( y ):
    global x,N
    store = np.zeros( x.size, dtype = complex)
    for k in range(-N,N+1):
        store += y[k+N]*np.exp(1j*k*x)
    return store

def get_FP_op(N):
    return get_ddx(N).dot( mult_by_sin(2) )

def get_H_op(N):
    FP_op = get_FP_op(N)
    return 0.5*FP_op - 0.5*FP_op.conj().transpose()

FP_op = get_FP_op(N)
H_op = get_H_op(N)

# Initialize our half-density and density to a uniform distribution
psi_0 = np.zeros(2*N+1,dtype=complex)
psi_0[N] = 1.0
rho_0 = np.zeros(2*N+1,dtype=complex)
rho_0[N] = 1.0

# Integrate solutions
n_frames = 100
t = np.linspace(0,1.5,n_frames)

#Solving for psi
psi_arr = np.zeros( ( n_frames , 2*N+1), dtype=complex )
integrator = ode( lambda t,x : H_op.dot(x) )
integrator.set_integrator('zvode',method='bdf')
integrator.set_initial_value( psi_0 , t[0] )
psi_arr[0] = psi_0
for k in range(1,n_frames):
    if integrator.successful():
        integrator.integrate(t[k])
        psi_arr[k] = integrator.y
    else:
        print "Integration unsucessful"

rho_arr = np.zeros( ( n_frames , 2*N+1) , dtype=complex)
integrator = ode( lambda t,x : FP_op.dot(x) )
integrator.set_integrator('zvode',method='bdf')
integrator.set_initial_value( rho_0 , t[0] )
rho_arr[0] = rho_0
for k in range(1,n_frames):
    if integrator.successful():
        integrator.integrate(t[k])
        rho_arr[k] = integrator.y
    else:
        print "Integration unsucessful"

Louiville_L1 = np.zeros(n_frames)
quantum_L1 = np.zeros(n_frames)

error_H = np.zeros(n_frames)
error_fp = np.zeros(n_frames)

for k in range(n_frames):
    #plt.plot(x, ift(rho).real, 'b', x, np.zeros(len(x)) , 'k' )
    rho_exact = (np.exp(2*t[k])*np.sin(x)**2 + np.exp(-2*t[k])*np.cos(x)**2)**(-1)
    rho = rho_arr[k]
    psi = psi_arr[k]
    rho_spatial = ift(rho)
    psi_spatial = ift(psi)
    Louiville_L1[k] = np.mean(abs(rho_spatial))
    quantum_L1[k] = (np.abs(psi)**2).sum()
    error_H[k] = np.mean( np.abs( np.abs(psi_spatial)**2 - rho_exact) )
    error_fp[k] = np.mean( abs( rho_spatial - rho_exact) )
    #U = expm( t[k]*H_op)
    #psi = ift(np.dot( U , psi_0))
    #plt.plot(x, abs(psi)**2 , 'b', x, np.zeros(len(x)) , 'k' )
    #plt.plot(x,rho_exact,'r')
    #plt.grid(True)
    #plt.xlabel('x')
    #plt.title('t={:f}'.format(t[k]))
    #plt.axis([-np.pi , np.pi , -1. , 10.])
    #fname = './figures/figure_{:d}.png'.format(k)
    #plt.savefig(fname)
    #plt.clf()

#plotting final condition
plot_final_distributions = False
if plot_final_distributions:
    psi_spatial = ift(psi_arr[-1])
    rho_spatial[0] = 0
    rho_spatial[-1] = 0
    psi_spatial[0] = 0
    psi_spatial[-1] = 0
    rho_exact[0] = 0
    rho_exact[-1] = 0

    plt.fill( x , rho_exact , color='0.5')
    #plt.xlabel('x'), plt.ylabel('mass density',fontsize=15)
    plt.axis([-3.14,3.14,-1,20])
    plt.grid()
    plt.tight_layout()
    plt.show()

    plt.fill( x , rho_spatial , color='0.5' )
    #plt.xlabel('x'), plt.ylabel('mass density',fontsize=15)
    plt.axis([-3.14,3.14,-1,20])
    plt.grid()
    plt.tight_layout()
    plt.show()

    plt.fill( x , abs(psi_spatial)**2 , color='0.5')
    #plt.xlabel('x'), plt.ylabel('mass density',fontsize=15)
    plt.axis([-3.14,3.14,-1,20])
    plt.grid()
    plt.tight_layout()
    plt.show()

    #plotting L1 norm
    plt.plot( t , Louiville_L1, 'k-')
    plt.plot( t , quantum_L1, 'k-.')
    plt.xlabel('time',fontsize=15)
    plt.ylabel('$L^1$ norm')
    plt.grid()
    plt.tight_layout()
    plt.show()
    quit()


#NOW WE MAKE A CONVERGENCE PLOT
make_convergence_plot = False

if make_convergence_plot:
    t_final = 1.0
    n_trials = 70
    error_H = np.zeros(n_trials)
    error_FP = np.zeros(n_trials)
    for k in range(n_trials):
        #Solving for psi
        power = k+2
        N = int(1.1**power)
        psi_0 = np.zeros( 2*N+1 , dtype = complex )
        psi_0[N] = 1.
        H_op = get_H_op(N)
        integrator = ode( lambda t,x : H_op.dot(x) )
        integrator.set_integrator('zvode',method='bdf')
        integrator.set_initial_value( psi_0 , 0 )
        integrator.integrate(t_final)
        psi = integrator.y
        
        #Solving for rho
        rho_0 = np.zeros( 2*N+1 , dtype = complex )
        rho_0[N] = 1.
        FP_op = get_FP_op(N)
        integrator = ode( lambda t,x : FP_op.dot(x) )
        integrator.set_integrator('zvode',method='bdf')
        integrator.set_initial_value( rho_0 , 0 )
        integrator.integrate(t_final)
        rho = integrator.y
        
        #Compute errors
        rho_exact = (np.exp(2*t_final)*np.sin(x)**2 + np.exp(-2*t_final)*np.cos(x)**2)**(-1)
        rho_spatial = ift(rho)
        psi_spatial = ift(psi)
        #plt.plot( x , rho_exact,'k')
        #plt.plot( x , rho_spatial.real,'b')
        #plt.plot( x , psi_spatial.real**2,'r')
        #plt.show()
        error_H[k] = 2*np.pi*np.mean( np.abs( np.abs(psi_spatial)**2 - rho_exact) )
        error_FP[k] = 2*np.pi*np.mean( np.abs( rho_spatial - rho_exact) )

    res = 1.1**(np.arange(n_trials)+2)
    print res
    print error_H
    error_wavelets = np.array( [2.82689808601994,
        1.12546699460215,
        0.1980847118894,
        0.0139577038454514,
        0.000650046703770487,
        0.000223111162371991,
        0.00014807650626467])
    res_wavelets = np.array([56,88,152,280,536,1048,2072])
    plt.loglog( res, error_H,'k-.' )
    plt.loglog( res, error_FP, 'k-' )
    #plt.loglog( res_wavelets, error_wavelets, 'k' )
    plt.grid(True)
    #plt.title('L^1 error')
    plt.tight_layout()
    plt.show()
    quit()

#Let's test conservation of scalar products!
#Let function_1 = sin(x)
#Let function_2 = cos(x)

FS_1 = np.zeros( 2*N+1, dtype = complex)
FS_1[1+N] = 1.0 / (2j)
FS_1[-1+N] = -1.0 / (2j)
FS_2 = np.zeros( 2*N+1, dtype = complex)
FS_2[1+N] = 1.0 / (2)
FS_2[-1+N] = 1.0 / (2)
FS_3 = np.zeros( 2*N+1, dtype = complex)
FS_3[2+N] = 0.5 / (2j)
FS_3[-2+N] = -0.5/(2j)

x_grid = np.linspace(0,1,2*N+1)

def get_convolution_matrix( FS ):
    #produces the convolution matrix associated to a given Fourier series
    global N
    store = FS.nonzero()
    data = np.zeros( [store[0].size ,2*N+1] , dtype=complex)
    offsets = np.zeros( store[0].size )
    D = 2*N+1
    i = 0
    for k in FS.nonzero()[0]:
        data[i,:] = FS[k]*np.ones(D)
        offsets[i] = -(k-N)
        i += 1
    return sparse.dia_matrix( (data,offsets), shape = (D,D),dtype=complex )

cnv1 = get_convolution_matrix( FS_1)
cnv2 = get_convolution_matrix( FS_2)
cnv3 = get_convolution_matrix( FS_3)
#D = 2*N+1
#sigma = 0.1
#data = np.zeros([1,D])
#data[0,:] = np.exp( - sigma**2 * (np.arange(-N,N+1))**2)
#offsets = np.zeros(1)
#Smooth = sparse.dia_matrix( (data,offsets), shape=(D,D),dtype=complex)
#data[0,:] = np.exp(  sigma**2 * (np.arange(-N,N+1))**2)
#Sharpen = sparse.dia_matrix( (data,offsets), shape=(D,D),dtype=complex)

#cnv0 = Smooth.dot( cnv1.dot( Sharpen) )
#cnv1 = cnv0.copy()

from scipy.sparse import linalg

#EVOLVE cnv1 and cnv2
def bracket( A , B ):
    return A.dot(B) - B.dot(A)
t=0
dt = 0.01
t_final = 1.5
plot_norms = True
if plot_norms:
    operator_norm = []
    #u,s,v = linalg.eigsh( cnv1 , k=1)
    operator_norm.append( np.linalg.norm( cnv1.todense(),ord=2) )
while t < t_final:
    k1 = bracket( H_op , cnv1)
    k2 = bracket( H_op , cnv1 + dt*k1 / 2)
    k3 = bracket( H_op , cnv1 + dt*k2 / 2)
    k4 = bracket( H_op , cnv1 + dt*k3 )
    cnv1 += dt*(k1+2*k2+2*k3+k4) / 6.
    #k1 = bracket( H_op , cnv2 )
    #k2 = bracket( H_op , cnv2 + dt*k1 / 2)
    #k3 = bracket( H_op , cnv2 + dt*k2 / 2)
    #k4 = bracket( H_op , cnv2 + dt*k3 )
    #cnv2 += dt*(k1+2*k2+2*k3+k4) / 6.
    #k1 = bracket( H_op , cnv3 )
    #k2 = bracket( H_op , cnv3 + dt*k1 / 2)
    #k3 = bracket( H_op , cnv3 + dt*k2 / 2)
    #k4 = bracket( H_op , cnv3 + dt*k3 )
    #cnv3 += dt*(k1+2*k2+2*k3+k4) / 6.
    if plot_norms:
        #u,s,v = linalg.svds( cnv1 , k=1)
        operator_norm.append( np.linalg.norm( cnv1.todense(),ord=2 ) )
    t += dt

#EVOLVE FS_1 and FS_2 and FS_3
dfdt = lambda f : -FP_op.transpose().conj().dot(f)
t=0
if plot_norms:
    sup_norm = []
    sup_norm.append( ift(FS_1).real.max() )
t_array = []
t_array.append(t)
while t < t_final:
    k1 = dfdt( FS_1 )
    k2 = dfdt( FS_1 + dt*k1 / 2)
    k3 = dfdt( FS_1 + dt*k2 / 2)
    k4 = dfdt( FS_1 + dt*k3 )
    FS_1 += dt*(k1+2*k2+2*k3+k4) / 6.
    k1 = dfdt( FS_2 )
    k2 = dfdt( FS_2 + dt*k1 / 2)
    k3 = dfdt( FS_2 + dt*k2 / 2)
    k4 = dfdt( FS_2 + dt*k3 )
    FS_2 += dt*(k1+2*k2+2*k3+k4) / 6.
    k1 = dfdt( FS_3 )
    k2 = dfdt( FS_3 + dt*k1 / 2)
    k3 = dfdt( FS_3 + dt*k2 / 2)
    k4 = dfdt( FS_3 + dt*k3 )
    FS_3 += dt*(k1+2*k2+2*k3+k4) / 6.
    if plot_norms:
        sup_norm.append( ift(FS_1).real.max() )
    t += dt
    t_array.append(t)

#PLOT norms
if plot_norms:
    plt.plot( t_array, sup_norm ,'k-')
    plt.plot( t_array, operator_norm,'k-.' )
    plt.xlabel('time')
    plt.ylabel('sup/operator norm')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

#PLOT cnv1,cnv2 and cnv1*cnv2
from scipy.linalg import eig
y,psi = eig( cnv1.todense() )
rho_array = np.zeros( [ x.size , y.size] )
for k in range(y.size):
    rho_array[:,k] = np.abs(ift( psi[:,k] ))**2
    rho_array[:,k] = rho_array[:,k] / rho_array[:,k].max()

#This is the exact solution to func 1
denominator = np.sqrt( np.cos(x)**2 + np.exp(4*t_final)*np.sin(x)**2)
exact1= np.exp( 2*t_final ) * np.sin(x) / denominator

#This is the exact solution to func 2
denominator = np.sqrt( np.sin(x)**2 + np.exp(-4*t_final)*np.cos(x)**2)
exact2= np.exp( -2*t_final ) * np.cos(x) / denominator

#This is the exact solution to func 3
denominator = np.sqrt( np.cos(2*x)**2 + np.exp(4*t_final)*np.sin(2*x)**2)
exact3 = exact1*exact2

segments = []
colors = []
for i in range(y.size):
    for j in range(x.size-1):
        segments.append( [ (x[j],y[i]),(x[j+1],y[i])] )
        c = 1-rho_array[j,i]
        colors.append( (c,c,c) )

from matplotlib.collections import LineCollection
lc = LineCollection( segments , cmap='Reds', colors=colors,
        linewidths=2)
f,ax = plt.subplots()
ax.add_collection(lc)
plt.axis( [-np.pi, np.pi , -1.5 , 1.5 ] )
plt.grid(True)
plt.show()

plt.plot( x , exact1,'k')
plt.axis( [-np.pi, np.pi , -1.5 , 1.5 ] )
plt.grid(True)
plt.show()


Koopman_am = ift(FS_2).real * ift(FS_1).real
Koopman_ma = ift(FS_3).real
discrepency = np.abs( Koopman_am - Koopman_ma)
plt.plot( x , discrepency , 'k-')
plt.axis([-np.pi,np.pi,0.,0.08])
plt.grid(True)
plt.show()

#PLOT FS_1,FS_2, and FS_1*FS_1 and FS_3
plt.plot( x , ift(FS_1).real , 'b')
plt.grid(True)
plt.axis( [-np.pi, np.pi , -1.5 , 1.5 ] )
plt.show()
