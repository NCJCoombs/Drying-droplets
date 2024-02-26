# Drying-droplets

## Problem description

This code simulates the evaporation of sessile particle laden droplets using a continuum model which is solved via the finite element method.

The user must first have an installation of the finite element library [oomph-lib](https://oomph-lib.github.io/oomph-lib/doc/html/).

Though the driver codes supplied are specific to droplets with elliptical and triangular contact lines, modifications for more general geometries can be easily made. The dimensional physical parameters for the problem, denoted by an asterisk, are:

Parameter  | Symbol
------------- | -------------
Characteristic length scale | $L^*$
Initial droplet volume | $V_0^*$
Fluid density  | $\rho^*$
Fluid viscosity  | $\mu^*$
Surface tension  | $\sigma^*$
Permeability of jammed particles | $k^*$
Evaporation flux density | $\mathcal{J}^*$
Solutal diffusion coefficient | $D^*$

The dimensionless qunatities to be specified in the code, which may be passed as command line arguments, are:

Dimensionless parameter | Symbol | Formula | Command line argument
------------- | ------------- | ------------- | -------------
Pecl&#233;t number | $\mathrm{Pe}$ |  $L^* \mathcal{J}^* /(\rho^* \epsilon D^*)$ | ```--peclet```
Scaled capillary number | $\widetilde{\mathrm{Ca}}$ | $\mu^* \mathcal{J}^* / (\rho^* \sigma^* \epsilon^4)$ | ```--capillary```
Scaled inverse pore size | $\nu$ | $\epsilon L^* / \sqrt{k^* \mu^*}$ | ```--nu```
Initial solute volume fraction | $\phi_0$ | - | ```--phi_initial```

where $\epsilon = (2/\pi)V_0^*/{L^*}^3$ is the drop's aspect ratio. For the elliptical drop driver code, the user must also specify:

Parameter  | Command line argument
------------- | -------------
Major axis length, $a$ | ```--major_axis```
Minor axis length, $b$ | ```--minor_axis```
Evaporation mode | ```--evaporation_flag```
Theoretical "dry-out" time, $t_f$ | ```--t_f```
Initial drop apex height | ```--initial_height```

The evaporation mode may be kinetic of diffusive, and can be specified via ```--evaporation_flag 0``` and ```--evaporation_flag 1``` respectively.

For computationally intensive runs, the simulation may be terminated and restarted at a later date via the command line flag ```--restart 1```. This reads in data from the file ```restart.dat``` which has been created during the previous run.

## Installation

The code has been written for and successfully compiled with [version 2.0.0](https://github.com/oomph-lib/oomph-lib/releases/tag/v2.0.0) of oomph-lib. The directories for the elements and driver codes mimic the directory structure of oomph-lib. Once these files have been placed in the relevant locations, the user may run ```./autogen.sh``` to generate the makefile for each driver.

For 2D simulations we recommend using the [mumps](https://mumps-solver.org/index.php) linear solver rather than oomph-lib's default, SuperLU.

## Example scripts and output

Example bash scripts have been provided for drops with an elliptical and triangular contact line. For the former, the fields $(h,p,\phi)$ are outputted along the major and minor axis at every time step. At a given time increment (0.025 in the code), a full output of all fields $(x,y,h,p,\phi,u,v)$ is done on both the computational mesh and a fine uniform mesh for aid with post-processing. At timestep ```i``` these are outputted to ```full_solni.dat``` and ```full_soln_uniformi.dat``` respectively.
