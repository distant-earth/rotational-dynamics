import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 20})
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

t = []
P = []
new_angle = []
with open('RESULT') as table:
	line1 = table.readline()
	ecc = float(line1.split()[0])
	q = float(line1.split()[1]) #earthRad
	phi0 = np.rad2deg(float(line1.split()[2])) # degrees
	P0 = float(line1.split()[3]) # hours
	AC = float(line1.split()[4])
	BC = float(line1.split()[5])
	for line in table:
		t.append(float(line.split()[0])) #hours
		P.append(float(line.split()[1])) #hours
		new_angle.append(np.rad2deg(float(line.split()[2]))) # degrees
t = np.array(t)
P = np.array(P)
new_angle = np.array(new_angle)
		
def fmt_y(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if (b!=0):
    	return r'${} \cdot 10^{{{}}}$'.format(a, b)
    else:
    	return r'${}$'.format(a, b)

if not os.path.exists('./pictures'):
		os.makedirs('./pictures')

plt.figure(figsize=(8,8))
plt.title(fr'$P_0 = {round(P0, 2)}$ h, $\quad\gamma_0 = {round(phi0)}$'+r'$^{\circ}$'+'\n'+fr'$q$ = {round(q, 2)} au, '+
fr'$\quad e$ = {round(ecc, 2)}, '+'\n'+fr'$A$/$C = {round(AC, 2)}$, $\quad B$/$C = {round(BC, 2)}$', fontsize=18)
plt.xlabel(r'$t$, h')
plt.ylabel(r'$\Delta P$, h')
plt.gca().xaxis.set_major_locator(MultipleLocator(15))
plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
plt.gca().yaxis.set_major_formatter(fmt_y)
step = 1
plt.plot(t[::step], P[::step] - P0, 'r-', lw = 1)
plt.tight_layout()
plt.savefig('pictures/dP.eps', format = 'eps')

plt.figure(figsize=(8,8))
plt.title(fr'$P_0 = {round(P0, 2)}$ h, $\quad\gamma_0 = {round(phi0)}$'+r'$^{\circ}$'+'\n'+fr'$q$ = {round(q, 2)} au, '+
fr'$\quad e$ = {round(ecc, 2)}, '+'\n'+fr'$A$/$C = {round(AC, 2)}$, $\quad B$/$C = {round(BC, 2)}$', fontsize=18)
plt.xlabel(r'$t$, h')
plt.ylabel(r'$\Delta \gamma, ^\circ$')
plt.gca().xaxis.set_major_locator(MultipleLocator(15))
plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
plt.gca().yaxis.set_major_formatter(fmt_y)
step = 1
plt.plot(t[::step], new_angle[::step] - phi0, 'r-', lw = 1)
plt.tight_layout()
plt.savefig('pictures/d_gamma.eps', format = 'eps')

plt.show()
