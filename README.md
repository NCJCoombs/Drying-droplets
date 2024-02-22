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

The dimensionless qunatities to be specified in the code are:

Dimensionless parameter | Symbol | Formula
------------- | ------------- | -------------
Pecl&#233;t number | $\mathrm{Pe}$ |  $L^* \mathcal{J}^* /(\rho^* \epsilon D^*)$
Scaled capillary number | $\widetilde{\mathrm{Ca}}$ | $\mu^* \mathcal{J}^* / (\rho^* \sigma^* \epsilon^4)$
Scaled inverse pore size | $\nu$ | $\epsilon L^* / \sqrt{k^* \mu^*}$

where $\epsilon = (2/\pi)V_0/L^3$ is the drop's aspect ratio.

