function [ddicdt,ddicds,ddicdalk,ddicdpco2]=calc_ddicdvar(grid,tin,sin,alkin,pin,co2in)
% [ddicdth,ddicds,ddicdalk,ddicdpco2]=calc_ddicdvar(grid,surf_theta,surf_salt,surf_alk.*molperm3_2_umolperkg,surf_po4.*molperm3_2_umolperkg,atm_pco2);
%
% Calculate the dDIC/d? coefficients from surface data
%
% Remember, CO2SYS uses umol/kg inputs

hfsurf=grid.hfacc(:,:,1);
tsurf=nanmean(tin,4).*hfsurf;
ssurf=nanmean(sin,4).*hfsurf;
psurf=nanmean(pin,4).*hfsurf;
asurf=nanmean(alkin,4).*hfsurf;

tmean=nansum(tsurf(:).*grid.rac(:))./nansum(grid.rac(:).*hfsurf(:));
smean=nansum(ssurf(:).*grid.rac(:))./nansum(grid.rac(:).*hfsurf(:));
amean=nansum(asurf(:).*grid.rac(:))./nansum(grid.rac(:).*hfsurf(:));
pmean=nansum(psurf(:).*grid.rac(:))./nansum(grid.rac(:).*hfsurf(:));

%simean=nansum(sisurf(:).*grid.rac(:))./nansum(grid.rac(:).*hfsurf(:));
% pco2mean=pco2;

% tmean=nansum(tsurf(:).*grid.volc(:).*grid.hfacc(:))./nansum(grid.volc(:).*grid.hfacc(:));
% smean=nansum(ssurf(:).*grid.volc(:).*grid.hfacc(:))./nansum(grid.volc(:).*grid.hfacc(:));
% amean=nansum(asurf(:).*grid.volc(:).*grid.hfacc(:))./nansum(grid.volc(:).*grid.hfacc(:));
% pmean=nansum(psurf(:).*grid.volc(:).*grid.hfacc(:))./nansum(grid.volc(:).*grid.hfacc(:));
pco2mean=nanmean(co2in(:));

pco2var=[1:50:3000]';
tvar=[floor(nanmin(tsurf(:))):1:ceil(nanmax(tsurf(:)))]';
svar=[floor(nanmin(ssurf(:))):0.2:ceil(nanmax(ssurf(:)))]';
avar=[floor(nanmin(asurf(:))):10:ceil(nanmax(asurf(:)))]';

% use climatological Si data for saturated carbon calculations,
% i think it only has a small effect.

if exist('~/Dropbox_Work/cfluxes_matlab/obs/woa13_annual_silicate.nc','file')
    % use local copy
    addr='~/Dropbox_Work/cfluxes_matlab/obs/woa13_annual_silicate.nc';
else
    % use remote copy
    addr='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA13/DATA/silicate/netcdf/all/1.00/woa13_all_i00_01.nc';
end

woa_si=nc_varget(addr,'i_an');  % umol/l
woa_latc=nc_varget(addr,'lat');
woa_lonc=nc_varget(addr,'lon'); %-179.5--179.5
woa_depth=nc_varget(addr,'depth');

% Rotate the atlantic centred data
%m_proj('mollweide','clong',0);
%[~,ri] = globe(woa_lonc,180);
ilon=[181:360,1:180];
[~,ilat,~]=intersect(woa_latc,grid.latc);
perl_2_perkg=1000./1024.5; % e.g. mol l-1 * 1000 l m-3 * 1/rho m3 kg-1 -> mol kg-1

[xin,yin,zin]=ndgrid(woa_lonc+180,woa_latc,woa_depth);
[xout,yout,zout]=ndgrid(grid.lonc,grid.latc,grid.zc);

Fi=griddedInterpolant(xin,yin,zin,inpaint_nans(permute(squeeze(woa_si(:,:,:,ilon)),[3,2,1]),4));
woa_si=Fi(xout,yout,zout).*grid.hfacc.*perl_2_perkg; % convert to umol/kg

sisurf=woa_si(:,:,1);
simean=nansum(sisurf(:).*grid.rac(:))./nansum(grid.rac(:));

%% See the end of this script for help deciphering the CO2SYS syntax.
% Calculate T dependence
A=CO2SYS(pco2mean,amean,4,1,smean,tvar,tvar,...
    0,0,simean,pmean,1,4,1);
tcsat=A(:,2); % umol/kg

[b,~,~,~,stats]=regress(tcsat,[ones(length(tcsat),1) tvar]);

figure
subplot(221)
scatter(tvar,tcsat,75,'filled')
hold on
plot([tvar(1) tvar(end)],[(b(1)+b(2)*tvar(1)) (b(1)+b(2)*tvar(end))],'r--','LineWidth',2)
ylabel('Csat (umol/kg)','FontSize',12)
xlabel('Temperature (degC)','FontSize',12)
title(['dCsat/dT = ',num2str(b(2)),' umol kg-1 K-1, R2 = ',num2str(stats(1))],'FontSize',14)
set(gca,'FontSize',12); box on; axis tight %orient landscape
%     print -dpsc2 ddicdt.ps
%     fixpslinestyle('ddicdt.ps')

ddicdt=b(2)*1e-6; % mol/kg/K

%% Calculate S dependence
A=CO2SYS(pco2mean,amean,4,1,svar,tmean,tmean,...
    0,0,simean,pmean,1,4,1);
scsat=A(:,2);

[b,~,~,~,stats]=regress(scsat,[ones(length(scsat),1) svar]);

subplot(222)
scatter(svar,scsat,75,'filled')
hold on
plot([svar(1) svar(end)],[(b(1)+b(2)*svar(1)) (b(1)+b(2)*svar(end))],'r--','LineWidth',2)
ylabel('Csat (umol/kg)','FontSize',12)
xlabel('Salinity','FontSize',12)
title(['dCsat/dS = ',num2str(b(2)),' umol kg-1 psu-1, R2 = ',num2str(stats(1))],'FontSize',14)
set(gca,'FontSize',12); box on; axis tight %orient landscape
%     print -dpsc2 ddicds.ps
%     fixpslinestyle('ddicds.ps')

ddicds=b(2)*1e-6; % mol/kg/psu

%% Calculate ALK dependence
A=CO2SYS(pco2mean,avar,4,1,smean,tmean,tmean,...
    0,0,simean,pmean,1,4,1);
acsat=A(:,2); % umol/kg

[b,~,~,~,stats]=regress(acsat,[ones(length(acsat),1) avar]);

subplot(223)
scatter(avar,acsat,75,'filled')
hold on
plot([avar(1) avar(end)],[(b(1)+b(2)*avar(1)) (b(1)+b(2)*avar(end))],'r--','LineWidth',2)
ylabel('Csat (umol/kg)','FontSize',12)
xlabel('Alkalinity (umol/kg)','FontSize',12)
title(['dCsat/dALK = ',num2str(b(2)),', R2 = ',num2str(stats(1))],'FontSize',14)
set(gca,'FontSize',12); box on; axis tight %orient landscape
%     print -dpsc2 ddicdalk.ps
%     fixpslinestyle('ddicdalk.ps')

ddicdalk=b(2); % umol/kg/(umol/kg)

%% Calculate pCO2 dependence
A=CO2SYS(pco2var,amean,4,1,smean,tmean,tmean,...
    0,0,simean,pmean,1,4,1);
pco2csat=A(:,2); % umol/kg
revelle=A(:,29);

%[b,~,~,~,stats]=regress(revelle,[ones(length(revelle),1) pco2var]);
ddicdpco2=(ones(size(pco2var)).*pco2csat)./(revelle.*pco2var);

%[b,~,~,~,stats]=regress(pco2csat(15:end),[ones(length(pco2csat(15:end)),1) pco2var(15:end)]);
rval=corr(pco2csat(2:end),(pco2csat(1:end-1)+ddicdpco2(2:end).*unique(diff(pco2var))));
subplot(224)
scatter(pco2var,pco2csat,75,'filled')
hold on
plot(pco2var(2:end),(pco2csat(1:end-1)+ddicdpco2(2:end).*unique(diff(pco2var))),'r--','LineWidth',2)
ylabel('Csat (umol/kg)','FontSize',12)
xlabel('pCO2 (uatm)','FontSize',12)
title({'dCsat/dpCO2 (umol kg-1 uatm-1) is variable,';['depending on pCO2. R2 = ',num2str(rval)]},'FontSize',14)
set(gca,'FontSize',12,'XLim',[0 3000]); box on;
t=text(750,1500,'$\frac{\delta C_{sat}}{\delta pCO_2}=\frac{\overline{C_{sat}}}{B pCO_2}$');
set(t,'Interpreter','Latex','FontSize',14);
orient landscape
print -dpsc2 ddicdvar.ps
fixpslinestyle('ddicdvar.ps')

ddicdpco2=[pco2var,ddicdpco2.*1e-6];
end

% help CO2SYS
%  **************************************************************************
%   CO2SYS is a MATLAB-version of the original CO2SYS for DOS.
%   CO2SYS calculates and returns the state of the carbonate system of 
%   oceanographic water samples, if supplied with enough input.
%   Also, it is useful for doing pH scale conversions. 
%  
%   Please note that his software is intended to be exactly identical to the 
%   DOS and Excel versions that have been released previously. Several coding
%   errors (typo's, etc.) were discovered in the initial versions of this 
%   matlab routine, which have all been corrected before version 1.00 was
%   released on CDIAC in June 2009. It is now considered to be fully conforming
%   to the version for Excel and DOS, meaning that results obtained should be
%   indentical for identical input.
%  
%   For lots more info please have a look at:
%   Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
%   CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
%   Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
%   Oak Ridge, Tennessee. 
%   http://cdiac.ornl.gov/oceans/co2rprt.html
%  **************************************************************************
%  
%   Syntax:
%   [DATA,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%         ...SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,...
%         ...K1K2CONSTANTS,KSO4CONSTANTS)
%  
%   Syntax example:
%   [a,headers,niceheaders]=CO2SYS(2400,2200,1,2,35,0,25,4200,0,15,1,1,4,1)
%  
%  **************************************************************************
%  
%   INPUT:
%  
%     PAR1  (some unit) : scalar or vector of size n
%     PAR2  (some unit) : scalar or vector of size n
%     PAR1TYPE       () : scalar or vector of size n (*)
%     PAR2TYPE       () : scalar or vector of size n (*)
%     SAL            () : scalar or vector of size n
%     TEMPIN  (degr. C) : scalar or vector of size n 
%     TEMPOUT (degr. C) : scalar or vector of size n 
%     PRESIN     (dbar) : scalar or vector of size n 
%     PRESOUT    (dbar) : scalar or vector of size n
%     SI    (umol/kgSW) : scalar or vector of size n
%     PO4   (umol/kgSW) : scalar or vector of size n
%     pHSCALEIN         : scalar or vector of size n (**)
%     K1K2CONSTANTS     : scalar or vector of size n (***)
%     KSO4CONSTANTS     : scalar or vector of size n (****)
%  
%    (*) Each element must be an integer, 
%        indicating that PAR1 (or PAR2) is of type: 
%    1 = Total Alkalinity
%    2 = DIC
%    3 = pH
%    4 = pCO2
%    5 = fCO2
%   
%    (**) Each element must be an integer, 
%         indicating that the pH-input (PAR1 or PAR2, if any) is at:
%    1 = Total scale
%    2 = Seawater scale
%    3 = Free scale
%    4 = NBS scale
%   
%    (***) Each element must be an integer, 
%          indicating the K1 K2 dissociation constants that are to be used:
%    1 = Roy
%    2 = Goyet & Poisson
%    3 = Hansson  refit by Dickson and Millero
%    4 = Mehrbach refit by Dickson and Millero (PREFERRED)
%    5 = Hansson and Mehrbach Refit by Dickson and Millero
%    6 = Geosecs
%    7 = Peng, Pure Water
%    8 = Millero, 1979 - DONT USE THIS - DOESN'T WORK
%    9 = Millero, 2006
%   10 = Lueker, 2000 
%    
%   
%    (****) Each element must be an integer, 
%           indicating the KSO4 dissociation constants that are to be used:
%    1 = Dickson's KSO4 (PREFERRED)
%    2 = Khoo's KSO4
%  
%  **************************************************************************%
%  
%   OUTPUT: * an array containing the following parameter values (one row per sample):
%           * a cell-array containing crudely formatted headers
%           * a cell-array containing nicely formatted headers
%  
%      POS  PARAMETER        UNIT
%  
%      01 - TAlk             (umol/kgSW)
%      02 - TCO2             (umol/kgSW)
%      03 - pHin             ()
%      04 - pCO2in           (uatm)
%      05 - fCO2in           (uatm)
%      06 - HCO3in           (umol/kgSW)
%      07 - CO3in            (umol/kgSW)
%      08 - CO2in            (umol/kgSW)
%      09 - BAlkin           (umol/kgSW)
%      10 - OHin             (umol/kgSW)
%      11 - PAlkin           (umol/kgSW)
%      12 - SiAlkin          (umol/kgSW)
%      13 - Hfreein          (umol/kgSW)
%      14 - RevelleFactorin  ()
%      15 - OmegaCain        ()
%      16 - OmegaArin        ()
%      17 - xCO2in           (ppm)
%      18 - pHout            ()
%      19 - pCO2out          (uatm)
%      20 - fCO2out          (uatm)
%      21 - HCO3out          (umol/kgSW)
%      22 - CO3out           (umol/kgSW)
%      23 - CO2out           (umol/kgSW)
%      24 - BAlkout          (umol/kgSW)
%      25 - OHout            (umol/kgSW)
%      26 - PAlkout          (umol/kgSW)
%      27 - SiAlkout         (umol/kgSW)
%      28 - Hfreeout         (umol/kgSW)
%      29 - RevelleFactorout ()
%      30 - OmegaCaout       ()
%      31 - OmegaArout       ()
%      32 - xCO2out          (ppm)
%      33 - pHin (Total)     ()          
%      34 - pHin (SWS)       ()          
%      35 - pHin (Free)      ()          
%      36 - pHin (NBS)       ()          
%      37 - pHout(Total)     ()          
%      38 - pHout(SWS)       ()          
%      39 - pHout(Free)      ()          
%      40 - pHout(NBS)       () 
%      41 - TEMPIN           ***    
%      42 - TEMPOUT          ***
%      43 - PRESIN           ***
%      44 - PRESOUT          ***
%      45 - PAR1TYPE         ***
%      46 - PAR2TYPE         ***
%      47 - K1K2CONSTANTS    ***
%      48 - KSO4CONSTANTS    *** 
%      49 - pHSCALEIN        ***
%      50 - SAL              ***
%      51 - PO4              ***
%      52 - SI               ***
%  
%      *** SIMPLY RESTATES THE INPUT BY USER 
%  
%  **************************************************************************
%  
%   CO2SYS originally by Lewis and Wallace 1998
%  
%   Converted to MATLAB by Denis Pierrot at
%   CIMAS, University of Miami, Miami, Florida
%   
%   Vectorization, internal refinements and speed improvements by
%   Steven van Heuven, University of Groningen, The Netherlands.
%  
%   Questions, bug reports et cetera: svheuven@gmail.com
%  
%  **************************************************************************
%  
%   This is version 1.01 (svhsvnrev62)
%   changes since 1.00:
%   - added a note explaining that all known bugs were removed before release of 1.00
%  
%   Uploaded to CDIAC ( http://cdiac.ornl.gov/oceans/co2rprt.html ) at June 11th, 2009.
%  
%  **************************************************************************