# Numerical integration of the rotational dynamics of an asteroid during close approach to the Earth

This repository provides the source code for studying the rotational dynamics of an asteroid during close approach to the Earth. The project is built upon the [DOP853 integrator](http://www.unige.ch/~hairer/prog/nonstiff/dop853.f) developed by E. Hairer and G. Wanner. See [this article](https://link.springer.com/article/10.1134/S0038094623700107) for further details on the problem.

## Folder overview:
By default, this folder contains the following files:
- **main.f90** — project main program
- **modd.f90** — Fortran module with functions and subroutines
- **dop853.f** — DOP853 integrator
- **plotting.py** — Python script allowing visualization of the results
- **Makefile**
- **INPUT** — file with a list of input parameters
- **LICENSE.md** and **README.md**

## INPUT:
You may specify the following parameters of the problem in the INPUT file:
- $P_0$ — initial rotational period (in hours)
- $q$ — pericentric distance of the hyperbolic orbit (in $R_{Earth}$)
- $e$ — eccentricity of the hyperbolic orbit ($e > 1$)
- $\varphi_0, \theta_0, \psi_0$ — Euler angles to define the initial orientation of the rotational axis of the asteroid (in degrees)
- $\dfrac{A}{C}, \dfrac{B}{C}$ — ratios of the asteroid's moments of inertia
- $R_\mathrm{sph}$ — numerical integration is performed within a geocentric sphere of this radius (in $R_{Earth}$)

## Makefile:
Makefile allows you to use the following commands (Linux terminal):
- `make res` to (re)compile and (re)run the project. Results are written into the RESULT file (generated automatically).
- `make plot` to visualize the results. Two images are generated, shown and saved into './pictures' folder (created automatically).
- `make clean` to remove temporary files, including RESULT file.

_Before using_ `make clean`, _make sure to save your RESULT file elsewhere (if necessary)!_




