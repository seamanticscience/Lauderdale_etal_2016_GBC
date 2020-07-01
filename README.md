# CO2 Flux diagnostic routines from Lauderdale_etal_2016_GBC
# 

<a href="https://doi.org/10.5281/zenodo.885500"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.885500.svg" alt="DOI"></a>
![GitHub release (latest by date)](https://img.shields.io/github/v/release/seamanticscience/Lauderdale_etal_2016_GBC?color=1b3370)
![GitHub last commit](https://img.shields.io/github/last-commit/seamanticscience/Lauderdale_etal_2016_GBC?color=f44323)
![GitHub License](https://img.shields.io/github/license/seamanticscience/Lauderdale_etal_2016_GBC?color=ffa500)
<a href="http://doi.org/10.1002/2016GB005400"><img src="http://img.shields.io/badge/paper-doi:%2F10.1002%2F2016GB005400-lightgrey?link=http://doi.org/10.1002/2016GB005400.svg" alt="Link to paper at http://doi.org/10.1002/2016GB005400"></a>

MATLAB routines for the CO2 flux diagnostics calculated using a
steady-state model from the paper "Quantifying the drivers of
ocean-atmosphere CO2 fluxes" by Jonathan Lauderdale, Stephanie
Dutkiewicz, Ric Williams and Mick Follows. Please cite this paper if you
find these tools of use in your research.

Run `cflux_diags.m` to process the supplied model output, in turn
`co2_flux_steadystate.m` is called to calculate the different components
and plot budgets and fluxes. The framework evaluates the interplay
between:
  1. Surface heat and freshwater fluxes that influence the
potential saturated carbon concentration, which depends on changes in
sea surface temperature, salinity and alkalinity, 
  1. A residual, disequilibrium flux influenced by upwelling and entrainment of
remineralized carbon- and nutrient-rich waters from the ocean interior,
as well as rapid subduction of surface waters, 
  1. Carbon uptake and export by biological activity as both soft tissue and carbonate, and 
  1. The effect on surface carbon concentrations due to freshwater
precipitation or evaporation. 

In a steady state simulation of a
coarse-resolution ocean circulation and biogeochemistry model, the sum
of the individually determined components is close to the known total
flux of the simulation.

There are a few options that can be changed in the `cflux_diags` driver
routine, use *mitgcm_fluxes* and set *vertmix* to blank to reproduce the
"online" fluxes (blue line in Figure 1d and the majority of the other
plots in our paper). Use *mitgcm_monthly* and set *vertmix* to
*diffusion* to calculate the "offline" flux (yellow line in Figure 1d).
Vary *kave* to change the depth of vertical integration.

The model fields in *model_kpp* are 12x monthly averages of a 3 degree
configuration of MITgcm - see our paper for more details.

There are two external dependencies: 
  1. `nc_isvar` for checking variables in netcdf files, which is part of the [`mexcdf/snctools` toolbox](http://mexcdf.sourceforge.net/)
  1. `CO2SYS` for solving the carbon system, which can be [obtained from CDIAC](https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/)

The MEX files, which enable interaction between *MATLAB* and *FORTRAN*
routines for speedy calculation, should probably be recompiled for your
system. They can then be called just like regular *MATLAB* functions (see
the examples below). Note that faults caused by these routines can cause
*MATLAB* to crash without notice. I tried to minimize this risk, but just
dont run them if you have thesis data that you have spent a long time
processing. They are compatible with 32- or 64-bit array sizes, which
defaults to different settings depending on your version of MATLAB
(groan). There are several checks that are run, including setting NaNs
to zero in masked regions. However, if you are experiencing problems,
you can compile the routines with the "-DDEBUG_MESSAGES" flag to get a
whole load of verbose output.

Any questions or comments, please get in contact!

# Installing MEX components
Check your mex installation and compiler ([see this Gist for help](https://gist.github.com/seamanticscience/b592fe74683f7e9dfc06914ca6536423))
```
>> mex -v -setup FORTRAN
```

- Offline advective flux code. 
```
Usage: [uflux,vflux,wflux]=mit_advect_flux(advscheme,timestep,dx,dy,dz,rax,ray,rac,mask,uvel,vvel,wvel,field);
>> mex advection_flux.F
```

- Offline explicit (GM) diffusive flux code.
```
Usage: [udiff,vdiff,wediff]=diffusion_flux(dx,dy,dz,umask,vmask,cmask,rax,ray,rac,kux,kuz,kvy,kvz,kwx,kwy,field);
>> mex diffusion_flux.F
```

- Offline implicit diffusive flux code.
``` 
Usage: widiff=implicit_diffusion_flux(timestep,iterations,rac,dzf,dzc,cmask,diffkz,field);
>> mex implicit_diffusion_flux.F 
```

- Routine to calculate vertically-integrated flux divergence.
```
Usage: [horz_div,vert_div]=calc_divergence(dz,vol,kave,mask,uflux,vflux,wflux);
>> mex calc_divergence.F
```

- Routine to calculate vertically-integrated flux divergence using a larger footprint.
```
Usage: [horz_div,vert_div]=calc_smooth_divergence(dz,vol,kave,mask,uflux,vflux,wflux);
>> mex calc_smooth_divergence.F
```

- Routine to integrate other properties vertically.
``` 
Usage: [integral_field]=integrate_vertically(dz,kave,mask,field);
>> mex integrate_vertically.F
```
