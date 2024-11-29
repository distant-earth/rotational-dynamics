import astropy.units as u
import numpy as np
from scipy.interpolate import CubicSpline as interp
from scipy.optimize import minimize_scalar as minimize
from scipy.optimize import curve_fit
import scipy.linalg
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.mplot3d import Axes3D
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 18})
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

if not os.path.exists('./pictures'):
		os.makedirs('./pictures')

#============================================
# Extracting ephemeris from file:

with open('horizons_results.txt') as table:
	for count, line in enumerate(table):
		if (line.strip() == '$$SOE'):
			N_start = count + 1
		if (line.strip() == '$$EOE'):
			N_end = count
			
time = 0
r = []
x = []
y = []
z = []

with open('horizons_results.txt') as table:
	data = table.readlines()
	for i in range(N_start + 2, N_end, 3):
		dist = (float(data[i][31:52]) * u.km).to(u.earthRad).value
		if (dist <= 100):
			r.append(dist)
			time += 1
			x.append((float(data[i-1][4:26]) * u.km).to(u.earthRad).value)
			y.append((float(data[i-1][30:52]) * u.km).to(u.earthRad).value)
			z.append((float(data[i-1][56:78]) * u.km).to(u.earthRad).value)

x = np.array(x)
y = np.array(y)
z = np.array(z)

#============================================
# Constructing best-fit orbit plane:

data = np.array([[x[i], y[i], z[i]] for i in range(len(x))])
mn = np.min(data, axis=0)
mx = np.max(data, axis=0)
# regular grid covering the domain of the data
X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
XX = X.flatten()
YY = Y.flatten()

# best-fit linear plane
A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
C,resid,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

# evaluate it on grid
Z = C[0]*X + C[1]*Y + C[2]

# points above and below the plane
dx = (z - C[1] * y - C[2]) / C[0] - x
dy = (z - C[0] * x - C[2]) / C[1] - y
dz = C[0] * x + C[1] * y + C[2] - z
Di = z - (C[0] * x + C[1] * y)
sgn = np.array(Di - C[2]) 
sgn = sgn / abs(sgn)
data_up = []
data_down = []
for i in range(len(data)):
	if (sgn[i] > 0):
		data_up.append(data[i,:])
	else:
		data_down.append(data[i,:])
data_up = np.array(data_up)
data_down = np.array(data_down)

# plot points and fitted surface
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(projection='3d')
ax.view_init(elev=27, azim=-135)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
ax.scatter(data_up[:,0], data_up[:,1], data_up[:,2], c='r', s=10)
ax.scatter(data_down[:,0], data_down[:,1], data_down[:,2], c='k', s=10)
ax.scatter(0, 0, 0, c='b', s=12)
ax.set_xlabel('\n'r'$X$, $R_E$')
ax.set_ylabel('\n'r'$Y$, $R_E$')
ax.set_zlabel('\n'r'$Z$, $R_E$')
ax.axis('auto')
ax.axis('tight')
plt.savefig('pictures/orbit.jpg', format = 'jpg', dpi = 600)

#============================================
# Constructing real orbit projection onto the plane:

z_fitted = C[0] * x + C[1] * y + C[2] # projection

t = np.arange(0, time, 1) # time in hours with 1-hour step
r_t = interp(t, r) # interpolation of r(t)

theta = []
theta.append(0.0)
for i in range(1, time):
	del_theta = np.arccos((x[i-1] * x[i] + y[i-1] * y[i] + z_fitted[i-1] * z_fitted[i]) / (np.linalg.norm((x[i-1],y[i-1],z_fitted[i-1])) * np.linalg.norm((x[i],y[i],z_fitted[i]))))
	theta.append(theta[i-1] + del_theta)

theta_t = interp(t, theta) # interpolation of theta(t)

tau, d = minimize(r_t).x, minimize(r_t).fun # pericentric time and distance

f = theta - theta_t(tau) # conversion from theta to f

def r_func(f, ecc):
	return d * (1 + ecc) / (1 + ecc * np.cos(f))

f_left = []
r_left = []
for i in range(len(f)):
	if f[i] <= 0.5: # not strictly 0 so that two hyperbolic branches intersect visibly on the plot
		f_left.append(f[i])
		r_left.append(r[i])
f_right = []
r_right = []
for i in range(len(f)):
	if f[i] >= -0.5:
		f_right.append(f[i])
		r_right.append(r[i])

# hyperbolic approximation
popt_l, pcov_l = curve_fit(r_func, xdata = f_left, ydata = r_left, p0=[4.0])
popt_r, pcov_r = curve_fit(r_func, xdata = f_right, ydata = r_right, p0=[4.0])
ecc_l = popt_l[0]
a_l = d / (ecc_l - 1)
ecc_r = popt_r[0]
a_r = d / (ecc_r - 1)
print('Before approach (f<0): \n a (R_E) = ', a_l, '\n e = ', ecc_l, '\n q (R_E) = ', a_l * (ecc_l - 1))
print('After approach (f>0): \n a (R_E) = ', a_r, '\n e = ', ecc_r, '\n q (R_E) = ', a_r * (ecc_r - 1))

# Plot distances between real orbit points and best-fit plane ('residuals')
plt.figure(figsize=(8,7))
plt.xlabel(r'$f$, rad')
plt.ylabel(r'$\Delta r, R_E$')
dx = (z - C[1] * y - C[2]) / C[0] - x
dy = (z - C[0] * x - C[2]) / C[1] - y
dz = C[0] * x + C[1] * y + C[2] - z
Di = z - (C[0] * x + C[1] * y)
sgn = np.array(Di - C[2]) 
sgn = sgn / abs(sgn)
plt.axhline(y = 0.0, lw = 0.8, color = 'grey', ls = '--')
plt.plot(f, sgn * np.sqrt(dx**2 + dy**2 + dz**2), 'r--', lw = 2)
plt.tight_layout()
plt.savefig('pictures/residuals.jpg', format = 'jpg', dpi=600)

# Plot hyperbolic approximation
plt.figure(figsize=(8,7))
plt.xlabel(r'$X$, $R_E$')
plt.ylabel(r'$Y$, $R_E$')
f = np.array(f)
r = np.array(r)
ff = f[r <= 100]
rr = r[r <= 100]
plt.plot(0, 0, marker="o", markersize=5, markeredgecolor="blue", markerfacecolor="blue")
plt.plot(rr * np.sin(ff), rr * np.cos(ff), 'k.', label = 'NASA JPL ephemeris')
plt.plot(r_func(f_left, popt_l[0]) * np.sin(f_left), r_func(f_left, popt_l[0]) * np.cos(f_left), 'r--', label = fr'Before approach ($e_1$ = {round(ecc_l, 2)})')
plt.plot(r_func(f_right, popt_r[0]) * np.sin(f_right), r_func(f_right, popt_r[0]) * np.cos(f_right), 'g--', label = fr'After approach ($e_2$ = {round(ecc_r, 2)})')
plt.legend(fontsize='12', loc='lower center', framealpha = 0)
plt.tight_layout()
plt.savefig('pictures/hyperb.jpg', format = 'jpg', dpi=600)

plt.show()
