import numpy as np
from qutip import *
import matplotlib.pyplot as plt

# 4 Feb 2017

#	-----------------------------------------------
#	parameters
gamma 	=	3								# nonlinearity parameter (large = quantum, small = classical)
kappa	=	5								# cavity damping
D 		=	1								# driving strength
g_r 	=	0.1								# coupling between vdP and atoms
g_c		=	0.1								# coupling of atoms to the cavity
N_vdp 	=	6								# dimension of the vdP
N_c		=	2								# dimension of the cavity

J1 = file_data_read('1_J_D1p00.dat')
J2 = file_data_read('2_J_D1p00.dat')
J3 = file_data_read('3_J_D1p00.dat')
J4 = file_data_read('4_J_D1p00.dat')
J5 = file_data_read('5_J_D1p00.dat')
J6 = file_data_read('6_J_D1p00.dat')
J7 = file_data_read('7_J_D1p00.dat')

Jind1 = file_data_read('1_Jind_D1p00.dat')
Jind2 = file_data_read('2_Jind_D1p00.dat')
Jind3 = file_data_read('3_Jind_D1p00.dat')
Jind4 = file_data_read('4_Jind_D1p00.dat')
Jind5 = file_data_read('5_Jind_D1p00.dat')
Jind6 = file_data_read('6_Jind_D1p00.dat')
Jind7 = file_data_read('7_Jind_D1p00.dat')

#	-----------------------------------------------
#	MAXIMA of < J_+J_- > and < J_+J_- >_{ind} 
#	-----------------------------------------------
N 			=	[1, 2, 3, 4, 5, 6, 7]
Jmax 		=	[0] * 7
Jmax[0]	=	max((J1[:,1].real).tolist())
Jmax[1]	=	max((J2[:,1].real).tolist())
Jmax[2]	=	max((J3[:,1].real).tolist())
Jmax[3]	=	max((J4[:,1].real).tolist())
Jmax[4]	=	max((J5[:,1].real).tolist())
Jmax[5]	=	max((J6[:,1].real).tolist())
Jmax[6]	=	max((J7[:,1].real).tolist())

Jindmax 		=	[0] * 7
Jindmax[0]	=	max((Jind1[:,1].real).tolist())
Jindmax[1]	=	max((Jind2[:,1].real).tolist())
Jindmax[2]	=	max((Jind3[:,1].real).tolist())
Jindmax[3]	=	max((Jind4[:,1].real).tolist())
Jindmax[4]	=	max((Jind5[:,1].real).tolist())
Jindmax[5]	=	max((Jind6[:,1].real).tolist())
Jindmax[6]	=	max((Jind7[:,1].real).tolist())


#	-----------------------------------------------
#	PLOT
#	-----------------------------------------------
fig = plt.figure(figsize=(17, 8))

plt.subplot(121)
plt.plot(J1[:,0], J1[:,1], lw='2',label='1 atoms')
plt.plot(J2[:,0], J2[:,1], lw='2',label='2 atoms')
plt.plot(J3[:,0], J3[:,1], lw='2',label='3 atoms')
plt.plot(J4[:,0], J4[:,1], lw='2',label='4 atoms')
plt.plot(J5[:,0], J5[:,1], lw='2',label='5 atoms')
plt.plot(J6[:,0], J6[:,1], lw='2',label='6 atoms')
plt.plot(J7[:,0], J7[:,1], lw='2',label='7 atoms')
plt.title(r'$N_{res}=%i$, $N_{cav}=%i$, $D=%.1f$, $g_r=%.2f$, $g_c=%.2f$, $\gamma_2/\gamma_1=%.1f$, $\kappa=%.1f$' % (N_vdp, N_c, D, g_r, g_c, gamma, kappa))
plt.xlabel(r'$\gamma_1 t$')
plt.ylabel(r'$\langle J_+J_-\rangle$')
plt.legend()

plt.subplot(122)
plt.plot(N, Jmax, 'o', label=r'$\langle J_+J_-\rangle^{max}$')
plt.plot(N, Jindmax, '^', label=r'$\langle J_+J_-\rangle_{ind}^{max}$')
plt.xlabel('Number of atoms')
plt.ylabel(r'$\langle J_+J_- \rangle^{max}$, $\langle J_+J_- \rangle_{ind}^{max}$')
plt.title(r'Maximum of $\langle J_+J_-\rangle$ and $\langle J_+J_-\rangle_{ind}$.')
plt.xlim([0.5,7.5])
plt.legend(loc='best')

plt.savefig('Jscaling_sync.pdf',Transparent=True, bbox_inches='tight')
plt.show()