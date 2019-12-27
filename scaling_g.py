import numpy as np
from qutip import *
import matplotlib.pyplot as plt

# 11 March 2017

#	-----------------------------------------------
#	parameters
Delta	=	0.1								# detuning
gamma 	=	3								# nonlinearity parameter (large = quantum, small = classical)
kappa	=	5								# cavity damping
D 		=	1								# driving strength
g_r 	=	0.1								# coupling between vdP and atoms
g_c		=	0.1								# coupling of atoms to the cavity
t 		=	np.linspace(0, 200, 1000)
N_vdp 	=	7								# dimension of the vdP
N_c		=	3								# dimension of the cavity


#	-----------------------------------------------
#	1 ATOM
#	-----------------------------------------------
n 		=	1
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(N_c))

Jm		=	sm1

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm])

file_data_store('1atom.dat', (np.vstack((t, result.expect[0]))).T)

print('1 atom')



#	-----------------------------------------------
#	2 ATOMS
#	-----------------------------------------------
n 		=	2
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(2), qeye(N_c))
sm2		=	tensor(qeye(N_vdp), qeye(2), sigmam(), qeye(N_c))

Jm		=	sm1 + sm2

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm])

file_data_store('2atoms.dat', (np.vstack((t, result.expect[0]))).T)

print('2 atoms')



#	-----------------------------------------------
#	3 ATOMS
#	-----------------------------------------------
n 		=	3
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(2,1), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(2), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(2), qeye(2), qeye(N_c))
sm2		=	tensor(qeye(N_vdp), qeye(2), sigmam(), qeye(2), qeye(N_c))
sm3		=	tensor(qeye(N_vdp), qeye(2), qeye(2), sigmam(), qeye(N_c))

Jm		=	sm1 + sm2 + sm3

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm])

file_data_store('3atoms.dat', (np.vstack((t, result.expect[0]))).T)

print('3 atoms')


#	-----------------------------------------------
#	4 ATOMS
#	-----------------------------------------------
n 		=	4
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm2		=	tensor(qeye(N_vdp), qeye(2), sigmam(), qeye(2), qeye(2), qeye(N_c))
sm3		=	tensor(qeye(N_vdp), qeye(2), qeye(2), sigmam(), qeye(2), qeye(N_c))
sm4		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), sigmam(), qeye(N_c))

Jm		=	sm1 + sm2 + sm3 + sm4

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm, Jind])

file_data_store('4atoms.dat', (np.vstack((t, result.expect[0]))).T)

print('4 atoms')



#	-----------------------------------------------
#	5 ATOMS
#	-----------------------------------------------
n 		=	5
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm2		=	tensor(qeye(N_vdp), qeye(2), sigmam(), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm3		=	tensor(qeye(N_vdp), qeye(2), qeye(2), sigmam(), qeye(2), qeye(2), qeye(N_c))
sm4		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), sigmam(), qeye(2), qeye(N_c))
sm5		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), sigmam(), qeye(N_c))

Jm		=	sm1 + sm2 + sm3 + sm4 + sm5

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm])

file_data_store('5atoms.dat', (np.vstack((t, result.expect[0]))).T)

print('5 atoms')


#	-----------------------------------------------
#	6 ATOMS
#	-----------------------------------------------
n 		=	6
g_new	=	g_c / np.sqrt(n)

psi0	=	tensor(basis(N_vdp,0), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(2,1), basis(N_c,0))

r 		=	tensor(destroy(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
c 		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), destroy(N_c))

sm1		=	tensor(qeye(N_vdp), sigmam(), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm2		=	tensor(qeye(N_vdp), qeye(2), sigmam(), qeye(2), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm3		=	tensor(qeye(N_vdp), qeye(2), qeye(2), sigmam(), qeye(2), qeye(2), qeye(2), qeye(N_c))
sm4		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), sigmam(), qeye(2), qeye(2), qeye(N_c))
sm5		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), sigmam(), qeye(2), qeye(N_c))
sm6		=	tensor(qeye(N_vdp), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), sigmam(), qeye(N_c))

Jm		=	sm1 + sm2 + sm3 + sm4 + sm5 + sm6

c_ops 	= 	[np.sqrt(gamma) * r ** 2, r.dag(), np.sqrt(kappa) * c]

H 		=	Delta * r.dag() * r + D * (r + r.dag()) + 1j * g_r * (r * Jm.dag() - r.dag() * Jm) + 1j * g_new * (c * Jm.dag() - c.dag() * Jm)
result 	= 	mesolve(H, psi0, t, c_ops, [Jm.dag() * Jm])

file_data_store('6atoms.dat', (np.vstack((t, result.expect[0]))).T)

print('6 atoms')