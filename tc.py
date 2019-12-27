import numpy as np
from qutip import *
import matplotlib.pyplot as plt


#	-----------------------------------------------
#	parameters
N 		=	10							# dimension of the oscillator
kappa   =   1							# damping of the qubits
g       =   0.1							# coupling of qubit to oscillator
Gamma	=	2 * g**2 / kappa			# single atom emission rate
omega_a =	1							# qubit frequency
omega_f =	1							# field frequency
tlist	=	np.linspace(0, 100, 1000)	# time


#	-----------------------------------------------
#	2 atoms
psi0_2	=	tensor(basis(N,0), basis(2,1), basis(2,1))		# initial state, field in vacuum and atoms excited

a_2 	=	tensor(destroy(N), qeye(2), qeye(2))
sz1_2	=	tensor(qeye(N), -sigmaz(), qeye(2))
sz2_2	=	tensor(qeye(N), qeye(2), -sigmaz())
sm1_2	=	tensor(qeye(N), destroy(2), qeye(2))
sm2_2	=	tensor(qeye(N), qeye(2), destroy(2))
Jm_2	=	sm1_2 + sm2_2
J_2_ind	=	sm1_2.dag() * sm1_2 + sm2_2.dag() * sm2_2		# independent atoms

H_2 		=	omega_f * a_2.dag() * a_2 + 0.5 * omega_a * (sz1_2 + sz2_2) + 1j * g * (a_2 * Jm_2.dag() - a_2.dag() * Jm_2)
c_ops_2 	= 	[np.sqrt(kappa) * a_2]
results_2	=	mesolve(H_2, psi0_2, tlist, c_ops_2, [Jm_2.dag() * Jm_2, J_2_ind])


#	-----------------------------------------------
#	3 atoms
psi0_3		=	tensor(basis(N,0), basis(2,1), basis(2,1), basis(2,1))		# initial state, field in vacuum and atoms excited

a_3 		=	tensor(destroy(N), qeye(2), qeye(2), qeye(2))
sz1_3		=	tensor(qeye(N), -sigmaz(), qeye(2), qeye(2))
sz2_3		=	tensor(qeye(N), qeye(2), -sigmaz(), qeye(2))
sz3_3		=	tensor(qeye(N), qeye(2), qeye(2), -sigmaz())
sm1_3		=	tensor(qeye(N), destroy(2), qeye(2), qeye(2))
sm2_3		=	tensor(qeye(N), qeye(2), destroy(2), qeye(2))
sm3_3		=	tensor(qeye(N), qeye(2), qeye(2), destroy(2))
Jm_3		=	sm1_3 + sm2_3 + sm3_3
J_3_ind		=	sm1_3.dag() * sm1_3 + sm2_3.dag() * sm2_3 + sm3_3.dag() * sm3_3		# independent atoms

H_3 		=	omega_f * a_3.dag() * a_3 + 0.5 * omega_a * (sz1_3 + sz2_3 + sz3_3) + 1j * g * (a_3 * Jm_3.dag() - a_3.dag() * Jm_3)
c_ops_3 	= 	[np.sqrt(kappa) * a_3]
results_3	=	mesolve(H_3, psi0_3, tlist, c_ops_3, [Jm_3.dag() * Jm_3, J_3_ind])


#	-----------------------------------------------
#	4 atoms
psi0_4		=	tensor(basis(N,0), basis(2,1), basis(2,1), basis(2,1), basis(2,1))		# initial state, field in vacuum and atoms excited

a_4 		=	tensor(destroy(N), qeye(2), qeye(2), qeye(2), qeye(2))
sz1_4		=	tensor(qeye(N), -sigmaz(), qeye(2), qeye(2), qeye(2))
sz2_4		=	tensor(qeye(N), qeye(2), -sigmaz(), qeye(2), qeye(2))
sz3_4		=	tensor(qeye(N), qeye(2), qeye(2), -sigmaz(), qeye(2))
sz4_4		=	tensor(qeye(N), qeye(2), qeye(2), qeye(2), -sigmaz())
sm1_4		=	tensor(qeye(N), destroy(2), qeye(2), qeye(2), qeye(2))
sm2_4		=	tensor(qeye(N), qeye(2), destroy(2), qeye(2), qeye(2))
sm3_4		=	tensor(qeye(N), qeye(2), qeye(2), destroy(2), qeye(2))
sm4_4		=	tensor(qeye(N), qeye(2), qeye(2), qeye(2), destroy(2))
Jm_4		=	sm1_4 + sm2_4 + sm3_4 + sm4_4
J_4_ind		=	sm1_4.dag() * sm1_4 + sm2_4.dag() * sm2_4 + sm3_4.dag() * sm3_4 + sm4_4.dag() * sm4_4	# independent atoms

H_4 		=	omega_f * a_4.dag() * a_4 + 0.5 * omega_a * (sz1_4 + sz2_4 + sz3_4 + sz4_4) + 1j * g * (a_4 * Jm_4.dag() - a_4.dag() * Jm_4)
c_ops_4 	= 	[np.sqrt(kappa) * a_4]
results_4	=	mesolve(H_4, psi0_4, tlist, c_ops_4, [Jm_4.dag() * Jm_4, J_4_ind])


#	-----------------------------------------------
#	plot
fig	=	plt.figure()

plt.plot(tlist, results_2.expect[0], lw='2', label='2 atoms')
plt.plot(tlist, results_2.expect[1], lw='2', ls='--', color='b')

plt.plot(tlist, results_3.expect[0], lw='2', label='3 atoms')
plt.plot(tlist, results_3.expect[1], lw='2', ls='--', color='g')

plt.plot(tlist, results_4.expect[0], lw='2', label='4 atoms')
plt.plot(tlist, results_4.expect[1], lw='2', ls='--', color='r')

plt.grid()
plt.legend()
plt.ylabel(r'$\langle J_{+}(t)J_{-}(t)\rangle$ (solid), $\langle J_{+}(t)J_{-}(t)\rangle_{ind}$ (dashed)')
plt.xlabel(r'$t$')
plt.title(r'Average squared polarization $\langle J_{+}(t)J_{-}(t)\rangle$.')
plt.show()

fig.savefig('superR.pdf',Transparent=True)