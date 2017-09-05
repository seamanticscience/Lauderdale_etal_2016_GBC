# CO2 Flux diagnostic routines from Lauderdale_etal_2016_GBC

MATLAB routines for the CO2 flux diagnostics calculated using a steady-state model from the paper "Quantifying the drivers of ocean-atmosphere CO2 fluxes" by Jonathan Lauderdale, Stephanie Dutkiewicz, Ric Williams and Mick Follows (http://dx.doi.org/10.1002/2016GB005400). Please cite this paper if you find these tools of use in your research.

Run cflux_diags.m to process the supplied model output, in turn co2_flux_steadystate.m is called to calculate the different components and plot budgets and fluxes. The framework evaluates the interplay between (1) surface heat and freshwater fluxes that influence the potential saturated carbon concentration, which depends on changes in sea surface temperature, salinity and alkalinity, (2) a residual, disequilibrium flux influenced by upwelling and entrainment of remineralized carbon- and nutrient-rich waters from the ocean interior, as well as rapid subduction of surface waters, (3) carbon uptake and export by biological activity as both soft tissue and carbonate, and (4) the effect on surface carbon concentrations due to freshwater precipitation or evaporation. In a steady state simulation of a coarse-resolution ocean circulation and biogeochemistry model, the sum of the individually determined components is close to the known total flux of the simulation. 

There are a few options that can be changed in the cflux_diags driver routine, use "mitgcm_fluxes" and set "vertmix" to blank to reproduce the "online" fluxes (blue line in Figure 1d and the majority of the other plots in our paper). Use "mitgcm_monthly" and set "vertmix" to "diffusion" to calculate the "offline" flux (yellow line in Figure 1d). Vary "kave" to change the depth of vertical integration.

The model fields in "model_kpp" are 12x monthly averages of a 3 degree configuration of MITgcm - see our paper for more details.

The MEX files, which enable interaction between MATLAB and Fortran routines for speedy calculation, should probably be recompiled for your system. They can then be called just like regular MATLAB functions (see the examples below). Note that faults caused by these routines can cause MATLAB to crash without notice. I tried to minimize this risk, but just dont run them if you have thesis data that you have spent a long time processing.

Any questions or comments, please get in contact!

https://zenodo.org/badge/98926569.svg

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
