MatLab routines for the CO2 flux diagnostics calculated using a steady-state model from the paper "Quantifying the drivers of ocean-atmosphere CO2 fluxes" by Jonathan Lauderdale, Stephanie Dutkiewicz, Ric Williams and Mick Follows (http://dx.doi.org/10.1002/2016GB005400)

The MEX files, which enable interaction between MatLab and Fortran routines should be recompiled for your system. Once done, they can be called just like regular MatLab functions

% Make object files from subroutines - these are modified snippets of MITgcm code
>> mex -c code/*.f

% Offline advective flux code.
>> mex mit_advect_flux.F gad_*.o
% Offline explicit (GM) diffusive flux code.
>> mex mit_diffuse_flux.F mit_diff_*.o
% Offline implicit diffusive flux code.
>> mex mit_difimpl_flux.F mit_implicit_rmix.o solve_tridiagonal.o
% Routine to calculate flux divergence.
>> mex mit_calc_divergence.F 
% Routine to integrate properties vertically.
>> mex mit_integrate_vert.F

