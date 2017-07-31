# CO2 Flux diagnostic routines from Lauderdale_etal_2016_GBC

MATLAB routines for the CO2 flux diagnostics calculated using a steady-state model from the paper "Quantifying the drivers of ocean-atmosphere CO2 fluxes" by Jonathan Lauderdale, Stephanie Dutkiewicz, Ric Williams and Mick Follows (http://dx.doi.org/10.1002/2016GB005400). Please cite this paper if you find these tools of use in your research.

The MEX files, which enable interaction between MATLAB and Fortran routines, should probably be recompiled for your system. They can then be called just like regular MATLAB functions (see the examples below).

Run cflux_diags.m to process the supplied model output and calls co2_flux_steadystate.m to calculate the different components and plot budgets and fluxes. There are a few options that can be changed in the cflux_diags driver routine, use "mitgcm_fluxes" and "vertmix" to blank to reproduce the "online" fluxes (blue line in Figure 1d and the majority of the other plots in our paper). Use "mitgcm_monthly" and set "vertmix" to "diffusion" to calculate the "offline" flux (yellow line in Figure 1d). Vary "kave" to change the depth of vertical integration.

The model fields in "model_kpp" are 12x monthly averages of a 3 degree configuration of MITgcm - see our paper for more details.

Any questions or comments, please get in contact!

# Installing MEX components
% Check your mex installation and compiler
>> mex -setup Fortran

% Make object files from subroutines - these are modified snippets of MITgcm code
>> mex -c code/*.f -outputdir code/

% Offline advective flux code.
% [uflux,vflux,wflux]=mit_advect_flux(timestep,advscheme,dx,dy,dz,rax,ray,rac,mask,uvel,vvel,wvel,scalar);
>> mex mit_advect_flux.F code/gad_*.o

% Offline explicit (GM) diffusive flux code.
% [udiff,vdiff,wediff]=mit_diffuse_flux(dx,dy,dz,umask,vmask,cmask,rax,ray,rac,kux,kuz,kvy,kvz,kwx,kwy,field);
>> mex mit_diffuse_flux.F code/mit_diff_*.o

% Offline implicit diffusive flux code.
% widiff=mit_difimpl_flux(timestep,iterations,rac,dzf,dzc,cmask,diffkz,field);
>> mex mit_difimpl_flux.F code/mit_implicit_rmix.o code/solve_tridiagonal.o
            
% Routine to calculate vertically-integrated flux divergence. 
% [horz_div,vert_div]=mit_calc_divergence(dz,vol,kave,mask,uflux,vflux,wtflux);
>> mex mit_calc_divergence.F 

% Routine to integrate properties vertically. 
% [integral_field]=mit_integrate_vert(dz,kave,mask,field);
>> mex mit_integrate_vert.F
