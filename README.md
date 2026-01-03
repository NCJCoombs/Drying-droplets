# Drying-droplets
Article: 
https://journals.aps.org/pre/abstract/10.1103/PhysRevE.110.064607

<img src="https://github.com/NCJCoombs/Drying-droplets/blob/main/birds_eye_elliptical_drops.png" width="800">

## Problem description

This code simulates the evaporation of sessile particle laden droplets using a continuum model which is solved via the finite element method.

The user must first have an installation of the finite element library [oomph-lib](https://oomph-lib.github.io/oomph-lib/doc/html/).

Though the driver codes supplied are specific to droplets with elliptical and triangular contact lines, modifications for more general geometries can be easily made. The governing equations may be solved on both structured and unstructured meshes (though structured meshes are recommended for efficiency). The dimensional physical parameters for the problem, denoted by an asterisk, are:

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

Dimensionless parameter | Symbol | Formula | Command line argument | Default value in code
------------- | ------------- | ------------- | -------------| -
Pecl&#233;t number | $\mathrm{Pe}$ |  $L^* \mathcal{J}^* /(\rho^* \epsilon D^*)$ | ```--peclet```| None (must be specified by user)
Scaled capillary number | $\widetilde{\mathrm{Ca}}$ | $\mu^* \mathcal{J}^* / (\rho^* \sigma^* \epsilon^4)$ | ```--capillary```| $10^{-5}$
Scaled inverse pore size | $\nu$ | $\epsilon L^* / \sqrt{k^* \mu^*}$ | ```--nu```|$1000$
Initial solute volume fraction | $\phi_0$ | - | ```--phi_initial```|$0.0256$

where $\epsilon = (2/\pi)V_0^* / {L^*}^3$ is the drop's aspect ratio. For the elliptical drop driver code, the user should also specify:

Parameter  | Command line argument | Default value in code
------------- | -------------| -
Major axis length, $a$ | ```--major_axis``` | $1.0$
Minor axis length, $b$ | ```--minor_axis```| $1.0$
Evaporation mode | ```--evaporation_flag``` | $0$
Theoretical "dry-out" time, $t_f$ | ```--t_f``` | None (must be specified by user)
Initial drop apex height | ```--initial_height```| None (must be specified by user)

The evaporation mode may be kinetic of diffusive, and can be specified via ```--evaporation_flag 0``` and ```--evaporation_flag 1``` respectively.

For computationally intensive runs, the simulation may be terminated and restarted at a later date via the command line flag ```--restart 1```. This reads in data from the file ```restart.dat``` which has been created during the previous run.

## File description
### ```user_src/thin_film_brinkman```
```(refineable_)thin_film_brinkman_elements.h```, ```(refineable_)thin_film_brinkman_elements.cc```: Calculates the bulk residual and Jacobian contributions for the non-linear system at the element level; see equations (A1-3) of the article. The refineable versions of these elements account for the presence of hanging nodes introduced by quadtree refinement.

```Tthin_film_brinkman_elements.h```, ```Tthin_film_brinkman_elements.cc```: Same as the above, but for triangular elements rather than quadrilateral.

### ```user_drivers/Thin_film_brinkman```
```thin_film_brinkman1D_refineable.cc```: Driver code specialized to axisymmetric drops. Used to obtain data for section III A 2 of the article.

```thin_film_brinkman2D_ellipse_structured.cc```: Driver code specialized to elliptical drops. Used to obtain data for section III B 2 of the article.

```thin_film_brinkman2D_triangle.cc```: Driver code specialized to a triangular drop with corners smoothed by a circular arc. Used to obtain data for section III C of the article.
## Installation

The code has been written for and successfully compiled with [version 2.0.0](https://github.com/oomph-lib/oomph-lib/releases/tag/v2.0.0) of oomph-lib. The directories for the elements and driver codes mimic the directory structure of oomph-lib. Once these files have been placed in the relevant locations, the user may run ```./autogen.sh``` to generate the makefile for each driver.

For 2D simulations we recommend using the [mumps](https://mumps-solver.org/index.php) linear solver rather than oomph-lib's default, SuperLU.

## Output and example run

At a given time increment (0.025 in the code), a full output of all fields $(x,y,h,p,\phi,u,v)$ is done. For the 2D driver codes, this is done on both the computational mesh and a fine uniform mesh for aid with post-processing. At timestep ```i``` these are outputted to ```full_solni.dat``` and ```full_soln_uniformi.dat``` respectively. The computational timesteps and simulation parameters are outputted to ```time.dat``` and ```parameters.dat``` respectively. for compatibility with post-processing software, the structure of the output may be altered by modifying the ```output``` function in ```thin_film_brinkman_elements.cc```.

As an example run, consider the following:

```./thin_film_brinkman1D_refineable --phi_c 0.02 --capillary 1e-5 --peclet 100 --evaporation_flag 1 --nu 1e3```

This command corresponds to an axisymmetric drop simulation with an initial solute volume fraction of $\phi_0 = 0.02$, $\widetilde{\text{Ca}}$ and $\nu$ values equal to their defaults, a Pecl&#233;t number of $\text{Pe} = 200$ and diffusive evaporation conditions.
