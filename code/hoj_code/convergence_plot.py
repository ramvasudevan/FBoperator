import numpy as np
from scipy.linalg import expm
#from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 32

ddx = np.zeros( [2*N+1,2*N+1] , dtype=complex)
for k in range(-N,N):
    ddx[k,k] = 1j*k

def mult_by_e( k ):
    return np.eye(2*N+1,k=-k,dtype=complex)

def mult_by_sin(k):
    return (mult_by_e(k) - mult_by_e(-k))/(2j)

def mult_by_cos(k):
    return (mult_by_e(k) + mult_by_e(-k))/2.

x = np.linspace(-np.pi, np.pi,200)
def ift( y ):
    global x
    store = 0j
    for k in range(-N,N):
        store += y[k]*np.exp(1j*k*x)
    return store

#X = sin(x) ddx

FP_op = np.dot( ddx , mult_by_sin(2) )
H_op = 0.5*FP_op - 0.5*FP_op.conj().transpose()

# Initialize our half-density and density to a uniform distribution
psi_0 = np.zeros(2*N+1,dtype=complex)
psi_0[0] = 1.0
rho_0 = psi_0**2

# Integrate solutions
n_frames = 90
t = np.linspace(0,2.0,n_frames)

Louiville_L1 = np.zeros(n_frames)
quantum_L1 = np.zeros(n_frames)

error_H = np.zeros(n_frames)
error_fp = np.zeros(n_frames)

for k in range(n_frames):
    #plt.plot(x, ift(rho).real, 'b', x, np.zeros(len(x)) , 'k' )
    rho_exact = (np.exp(2*t[k])*np.sin(x)**2 + np.exp(-2*t[k])*np.cos(x)**2)**(-1)
    rho = np.dot(expm( t[k] * FP_op ) , rho_0 )
    psi = np.dot(expm( t[k] * H_op ) , psi_0 )
    rho_spatial = ift(rho).real
    psi_spatial = ift(psi)
    Louiville_L1[k] = np.mean(abs(rho_spatial))
    quantum_L1[k] = (np.abs(psi)**2).sum()
    error_H[k] = np.mean( abs( np.abs(psi_spatial)**2 - rho_exact) )
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
psi_spatial = ift(psi)
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

