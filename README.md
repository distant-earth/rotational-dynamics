# Numerical integration of the rotational dynamics of an asteroid during close approach to the Earth

[![DOI](https://zenodo.org/badge/889539887.svg)](https://doi.org/10.5281/zenodo.14173662)

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
- **orbit** (subfolder)
  - **parser.py** — orbit approximation tool (new in version 2.0)

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

## Orbit approximation tool (new in version 2.0)

The project contains Python script that performs hyperbolic approximation of asteroid's orbit in the vicinity of the Earth (within a sphere of 100 Earth radii).

To use the script, do the following steps:
1) Generate ephemeris for the preferred asteroid using [Horizons System](https://ssd.jpl.nasa.gov/horizons/app.html#/).
Use the following settings:
- **Ephemeris Type:** vector table
- **Target Body:** your preferred asteroid
- **Coordinate Center:** Geocentric
- **Time Specification:** Specify time span. _We recommend a 2-4 days time span. For example, if close approach of an asteroid to the Earth is happening on the 2029-Apr-13 21:46 (Apophis, see [NASA small-body database](https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/)), consider using 2029-Apr-11 21:46 as start time and 2029-Apr-15 21:46 as end time._ 
Step size: 1 hour.
![изображение](https://github.com/user-attachments/assets/6d1e7ea9-5352-4105-9d70-682b87f04151)
- **Vector Table Settings:**
![изображение](https://github.com/user-attachments/assets/71c0f36a-7315-43a0-8f67-571591128c6b)
2) Save generated ephemeris as _'horizons_results.txt'_ (default name) into the **orbit** subfolder of this project.
3) Open the **orbit** subfolder.
4) From Linux terminal, run the following command: ```python3 parser.py```

The script generates three images, shows them and saves into 'orbit/pictures' folder (created automatically). Hyperbolic orbit parameters $a$, $e$ and $q = a(e - 1)$ are approximated separately for two hyperbolic branches and listed in the terminal output.
