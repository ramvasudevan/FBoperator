import numpy as np
from scipy import sparse
from scipy.integrate import ode
import matplotlib.pyplot as plt

N = 64

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

x = np.linspace(-np.pi, np.pi,1000)
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
psi_spatial = ift(psi_arr[-1])
rho_spatial[0] = 0
rho_spatial[-1] = 0
psi_spatial[0] = 0
psi_spatial[-1] = 0
rho_exact[0] = 0
rho_exact[-1] = 0

plt.fill( x , rho_exact )
plt.title('Exact'), plt.xlabel('x'), plt.ylabel('mass density')
plt.show()

plt.fill( x , rho_spatial )
plt.title('Standard Spectral'), plt.xlabel('x'), plt.ylabel('mass density')
plt.show()

plt.fill( x , abs(psi_spatial)**2 )
plt.title('GN Spectral'), plt.xlabel('x'), plt.ylabel('mass density')
plt.show()

#plotting L1 norm
plt.plot( t , Louiville_L1, 'b-')
plt.plot( t , quantum_L1, 'r-')
plt.xlabel('time')
plt.ylabel('$L^1$ norm')
plt.grid()
plt.show()


#NOW WE MAKE A CONVERGENCE PLOT
t_final = 1.0
n_trials = 15
error_H = np.zeros(n_trials)
error_FP = np.zeros(n_trials)
for k in range(n_trials):
    #Solving for psi
    power = k+2
    N = int(1.5**power)
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

res = 2**(np.arange(n_trials)+2)
print res
print error_H
plt.loglog( res, error_H,'r' )
plt.loglog( res, error_FP, 'b' )
plt.grid(True)
plt.title('L^1 error')
plt.show()


