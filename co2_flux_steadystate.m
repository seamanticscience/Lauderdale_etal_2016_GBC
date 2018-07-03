%
% Not a function, just a collection of decomposition decisions that I
% needed to make sure were consistent between the model and obs (ie I only 
% wanted to edit it once!)
% JML 11/11/14
%
% DIC BUDGET METHOD
% FCO2 = - div_adv_dic    % divergence of advection of carbon
%        + div_diff_dic   % divergence of diffusion of carbon 
%        + dic_bio_flux   % soft tissue
%        + dic_carb_flux  % carbonate
%        + Fw/rho.Cmean;  % freshwater fluxes
%
%STEADY STATE USING SURFACE FLUXES (as in Lauderdale et al 2016, GBC, doi:10.1002/2016GB005400.)
% FCO2 = gamm_T.Fheat/rho.Cp                            % heat fluxes
%      + Fw/rho.(gamm_S.Smean+gamm_A.Amean-Cmean)       % freshwater fluxes
%      - Rcp.(-div_adv_po4+div_diff_po4)h               % soft tissue
%      - 0.5.Rcaco3.Rcp.(-div_adv_po4+div_diff_po4)h    % carbonate
%      + (-div_adv_Cres+div_diff_Cres)h                 % disequilibrium
%
%% First, plot some budget terms
%% Plot budgets of temperature, salt, dic, alk and po4
landmask=change(grid.hfacc(:,:,1)','>',0.999,5);
landmask=change(landmask,'==',NaN,1);
landmask=change(landmask,'==',5,NaN);

plot_landmask='false'; b=[];

% Reduce the advective terms by order of magnitude to plot on same scale
fac=0.01;

% Temperature budget
figure
h(1)=subplot(331);
%[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(theta_adv_horz,3)'.*grid.hfacc(:,:,1)',-1e-4:1e-5:1e-4);canom;colormap(bluewhitered(20));caxis([-1e-4 1e-4]);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(theta_adv_horz,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Horizontal Advection of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
%[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(theta_adv_vert,3)'.*grid.hfacc(:,:,1)',-1e-4:1e-5:1e-4);canom;colormap(bluewhitered(20));caxis([-1e-4 1e-4]);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(theta_adv_vert,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Vertical Advection of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean((theta_adv_horz+theta_adv_vert),3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((theta_adv_horz+theta_adv_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Sum of Advection';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(theta_diff_horz,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_diff_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Horizontal Diffusion of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(theta_diff_vert,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Vertical Diffusion of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean((theta_diff_horz+theta_diff_vert),3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((theta_diff_horz+theta_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Sum of Diffusion';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(theta_surf_flux,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_surf_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Surface Flux of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(dthetadt,3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dthetadt,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Tendency of T';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean((theta_surf_flux+theta_adv_horz+theta_adv_vert+theta_diff_horz+theta_diff_vert),3)'.*grid.hfacc(:,:,1)',[-5e-4,-5e-5:5e-6:5e-5,5e-4]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((theta_surf_flux+theta_adv_horz+theta_adv_vert+theta_diff_horz+theta_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-5); axestop; end
title({'Sum of T components';'[K.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none'); 
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_budget_theta.ps'])

% Salt Budget
figure
h(1)=subplot(331);
%[~,c(1)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_horz,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_horz.*fac,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({['Horizontal Advection of S/',num2str(1/fac)];'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
%[~,c(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_vert,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_vert.*fac,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({['Vertical Advection of S/',num2str(1/fac)];'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean((salt_adv_horz+salt_adv_vert),3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_adv_horz+salt_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Sum of Advection';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(salt_diff_horz,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_diff_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Horizontal Diffusion of S';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(salt_diff_vert,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Vertical Diffusion of S';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean((salt_diff_horz+salt_diff_vert),3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((salt_diff_horz+salt_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Sum of Diffusion';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(salt_surf_flux,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_surf_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Surface Flux of S';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(dsaltdt,3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dsaltdt,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Tendency of S';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean((salt_surf_flux+salt_adv_horz+salt_adv_vert+salt_diff_horz+salt_diff_vert),3)'.*grid.hfacc(:,:,1)',[-5e-5,-5e-6:5e-7:5e-6,5e-5]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((salt_surf_flux+salt_adv_horz+salt_adv_vert+salt_diff_horz+salt_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-6); axestop; end
title({'Sum of S components';'[Psu.m/s]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_budget_salt.ps'])

%% Carbon Budget
figure
h(1)=subplot(331);
%[~,c(1)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_horz,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(dic_adv_horz.*fac,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({['Horizontal Advection of DIC/',num2str(1/fac)];'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
%[~,c(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_vert,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(dic_adv_vert.*fac,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({['Vertical Advection of DIC/',num2str(1/fac)];'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean((dic_adv_horz+dic_adv_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((dic_adv_horz+dic_adv_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of Advection';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(dic_diff_horz+dic_diff_vert,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_diff_horz+dic_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Horizontal+Vertical Diffusion of DIC';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(dic_co2_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_co2_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Surface Flux of CO2';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean((dic_carb_flux+dic_bio_flux),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((dic_carb_flux+dic_bio_flux),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Net Biogenic Flux';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Surface Dilution of DIC';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(ddicdt,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(ddicdt,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Tendency of DIC';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean((dic_vrt_flux+dic_co2_flux+dic_carb_flux+dic_bio_flux+dic_adv_horz+dic_adv_vert+dic_diff_horz+dic_diff_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((dic_vrt_flux+dic_co2_flux+dic_carb_flux+dic_bio_flux+dic_adv_horz+dic_adv_vert+dic_diff_horz+dic_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of DIC components';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_budget_dic.ps'])

% Alkalinity Budget
figure
h(1)=subplot(331);
%[~,c(1)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_horz,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(alk_adv_horz.*fac,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({['Horizontal Advection of ALK/',num2str(1/fac)];'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
%[~,c(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_vert,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(alk_adv_vert.*fac,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({['Vertical Advection of ALK/',num2str(1/fac)];'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean((alk_adv_horz+alk_adv_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((alk_adv_horz+alk_adv_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of Advection';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(alk_diff_horz,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_diff_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Horizontal Diffusion of ALK';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(alk_diff_vert,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Vertical Diffusion of ALK';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean((alk_bio_flux+alk_carb_flux),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((alk_bio_flux+alk_carb_flux),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Biological+Carbonate Flux';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(alk_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Surface Dilution of ALK';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(dalkdt,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dalkdt,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Tendency of ALK';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean((alk_vrt_flux+alk_carb_flux+alk_bio_flux+alk_adv_horz+alk_adv_vert+alk_diff_horz+alk_diff_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((alk_vrt_flux+alk_carb_flux+alk_bio_flux+alk_adv_horz+alk_adv_vert+alk_diff_horz+alk_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of ALK components';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_budget_alk.ps'])

%% Phosphate Budget
figure
h(1)=subplot(331);
%[~,c(1)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_horz,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(po4_adv_horz,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(po4_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Horizontal Advection of PO4';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
%[~,c(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_adv_vert,3)'.*grid.hfacc(:,:,1)',-5e-4:5e-5:5e-4);canom;colormap(bluewhitered(20));caxis([-5e-4 5e-4]);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(po4_adv_vert,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(po4_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Vertical Advection of PO4';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean((po4_adv_horz+po4_adv_vert),3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((po4_adv_horz+po4_adv_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of Advection';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(po4_diff_horz,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(po4_diff_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Horizontal Diffusion of PO4';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(po4_diff_vert,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(po4_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Vertical Diffusion of PO4';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean((po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of Diffusion';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(po4_bio_flux,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(po4_bio_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Biological Flux';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(dpo4dt,3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dpo4dt,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Tendency of PO4';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean((po4_bio_flux+po4_adv_horz+po4_adv_vert+po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[-2e-8,-2e-9:2e-10:2e-9,2e-8]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((po4_bio_flux+po4_adv_horz+po4_adv_vert+po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-5e-7); axestop; end
title({'Sum of PO4 components';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_budget_po4.ps'])
clear a b h

% Bio production
% if isempty(strfind(obsid,'ecco2darwin'));
%     biop=-(po4_adv_horz+po4_adv_vert+po4_diff_vert+po4_diff_horz);
% else
    % Add two important time varying terms
if exist('dpo4dt','var') 
    biop=-(po4_adv_horz+po4_adv_vert+po4_diff_vert+po4_diff_horz-dpo4dt);
else
    biop=-(po4_adv_horz+po4_adv_vert+po4_diff_vert+po4_diff_horz);
end

%% These decompositions are in the spirit of the original proposal, which neglected mixing
ugradtheta=-(theta_diff_horz+theta_diff_vert+theta_surf_flux);
ugradsalt =-( salt_diff_horz+ salt_diff_vert+ salt_surf_flux);

% if strcmpi(obsid,'mitgcm')
%     % use model fluxes of carbonate and bio
%     ugradalk=-(alk_vrt_flux+alk_carb_flux+alk_bio_flux+alk_diff_horz+alk_diff_vert);
% else
    ugradalk=-(alk_vrt_flux-(Rnp.*biop)+(Rcaco3.*Rcp.*biop)+alk_diff_horz+alk_diff_vert);
%end

ugradcsat=ugradtheta.*ddicdth+ugradsalt.*ddicds+ugradalk.*ddicdalk;

% It is still better to find the residual rather than use cres (maybe dpco2ddic missing?)
ugradcres=(dic_adv_horz+dic_adv_vert)-ugradcsat;

%% Alternatively, we probably know advection better than diffusion, so substitute those instead
gradkgradtheta=-(theta_adv_horz+theta_adv_vert+theta_surf_flux);
gradkgradsalt =-(salt_adv_horz +salt_adv_vert +salt_surf_flux );
gradkgradalk  =-( alk_vrt_flux-(Rnp.*biop)+(Rcaco3.*Rcp.*biop)+alk_adv_horz+alk_adv_vert);

gradkgradcsat =gradkgradtheta.*ddicdth+gradkgradsalt.*ddicds+gradkgradalk.*ddicdalk;
gradkgradcres =(dic_diff_horz+dic_diff_vert)-gradkgradcsat;

%% Next up, use calculated advection and diffusion for each term
advtheta=(theta_adv_horz+theta_adv_vert);
advsalt=(salt_adv_horz+salt_adv_vert);
advalk=(alk_adv_horz+alk_adv_vert);

csat_adv_horz=theta_adv_horz.*ddicdth+salt_adv_horz.*ddicds+alk_adv_horz.*ddicdalk;
csat_adv_vert=theta_adv_vert.*ddicdth+salt_adv_vert.*ddicds+alk_adv_vert.*ddicdalk;

cres_adv_horz=dic_adv_horz-csat_adv_horz;
cres_adv_vert=dic_adv_vert-csat_adv_vert;

advcsat=csat_adv_horz+csat_adv_vert; %advtheta.*ddicdth+advsalt.*ddicds+advalk.*ddicdalk;
advcres=cres_adv_horz+cres_adv_vert; %(dic_adv_horz+dic_adv_vert)-advcsat;

difftheta=(theta_diff_horz+theta_diff_vert);
diffsalt=(salt_diff_horz+salt_diff_vert);
diffalk=(alk_diff_horz+alk_diff_vert);

csat_diff_horz=theta_diff_horz.*ddicdth+salt_diff_horz.*ddicds+alk_diff_horz.*ddicdalk;
csat_diff_vert=theta_diff_vert.*ddicdth+salt_diff_vert.*ddicds+alk_diff_vert.*ddicdalk;

cres_diff_horz=dic_diff_horz-csat_diff_horz;
cres_diff_vert=dic_diff_vert-csat_diff_vert;

diffcsat=csat_diff_horz+csat_diff_vert; %difftheta.*ddicdth+diffsalt.*ddicds+diffalk.*ddicdalk;
diffcres=cres_diff_horz+cres_diff_vert; %(dic_diff_horz+dic_diff_vert)-diffcsat;

%% Cancel out advection and diffusion substitutions leaving just surface fluxes.
csat_surface_fluxes=-(theta_surf_flux.*ddicdth)-(salt_surf_flux.*ddicds)-(alk_vrt_flux.*ddicdalk); % Assuming preformed alkalinity here.
cres_surface_fluxes=(dic_adv_horz+dic_adv_vert+dic_diff_horz+dic_diff_vert)-csat_surface_fluxes;

%% Tendency of DIC is accounted for in the error in CRES, so subtract it
if exist('ddicdt','var')
    cres_adv_vert=cres_adv_vert-ddicdt;
    cres_surface_fluxes=cres_surface_fluxes-ddicdt;
end

%% Adjust the disequilibrium for outgassing of riverine DIC
 if strcmpi(obsid(1:4),'clim')
     tmp=grid.hfacc(:,:,1);
     river_flux=repmat(0.5e15/(12*60*60*24*365.25*nansum(grid.rac(:).*tmp(:))),size(dic_co2_flux));
     cres_adv_vert=cres_adv_vert+river_flux;
     cres_surface_fluxes=cres_surface_fluxes+river_flux;
     advcres=advcres+river_flux;
     ugradcres=ugradcres+river_flux;
     clear tmp river_flux
 end

%% Now recreate co2 fluxes from components...exchange different pools for actual fields

% The complete DIC budget
flux_co2_budget=dic_adv_horz+dic_adv_vert+dic_diff_horz+dic_diff_vert...
    +dic_carb_flux+dic_bio_flux+dic_vrt_flux;

% All the components in the spirit of original proposal (finding advection terms)
flux_co2_ugradcsat=ugradcsat+ugradcres+dic_diff_vert+dic_diff_horz... % substituted for dic_adv_horz+dic_adv_vert
    +(Rcp.*biop)+(Rcaco3.*Rcp.*biop)+dic_vrt_flux;

% All the components in the spirit of original proposal (finding diffusion terms)
flux_co2_gradkgradcsat=dic_adv_horz+dic_adv_vert+gradkgradcsat+gradkgradcres... % substituted for dic_diff_horz+dic_diff_vert
    +(Rcp.*biop)+(Rcaco3.*Rcp.*biop)+dic_vrt_flux;

% Substitution of both advection and diffusion terms
flux_co2_adv_diff=advcsat+advcres+diffcsat+diffcres... % substituted for dic_adv_horz+dic_adv_vert+dic_diff_horz+dic_diff_vert
    +(Rcp.*biop)+(Rcaco3.*Rcp.*biop)+dic_vrt_flux;
    
% Substitution using surface fluxes
flux_surface_fluxes=csat_surface_fluxes+cres_surface_fluxes...
    +(Rcp.*biop)+(Rcaco3.*Rcp.*biop)+dic_vrt_flux;

b=[];

%% Plot the total components compared to reference co2 fluxes
figure
h(1)=subplot(231);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(dic_co2_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(dic_co2_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'Reference CO2 fluxes';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(232);
[~,a(2)]=contourf(grid.lonc,grid.latc,-nanmean(flux_co2_budget,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,-nanmean(flux_co2_budget,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 fluxes from';'DIC budget [mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(233);
[~,a(3)]=contourf(grid.lonc,grid.latc,-nanmean(flux_co2_ugradcsat,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,-nanmean(flux_co2_ugradcsat,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 fluxes from';'ugradcsat [mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(234);
[~,a(4)]=contourf(grid.lonc,grid.latc,-nanmean(flux_co2_gradkgradcsat,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,-nanmean(flux_co2_gradkgradcsat,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 fluxes from';'gradkgradcsat [mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(235);
[~,a(5)]=contourf(grid.lonc,grid.latc,-nanmean(flux_co2_adv_diff,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,-nanmean(flux_co2_adv_diff,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 fluxes from';'adv and diff csat [mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(236);
[~,a(6)]=contourf(grid.lonc,grid.latc,-nanmean(flux_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on; if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,-nanmean(flux_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 fluxes from';'surface fluxes [mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_co2_flux_components.ps'])

%%
figure
ploth(1)=plot(grid.latc,squeeze(nansum(nanmean(dic_co2_flux,3).*grid.hfacc(:,:,1).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',2);
hold on
ploth(2)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_co2_budget,3).*grid.hfacc(:,:,1).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'g--','LineWidth',1);
ploth(3)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_co2_ugradcsat,3).*grid.hfacc(:,:,1).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'r','LineWidth',2);
ploth(4)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_co2_gradkgradcsat,3).*grid.hfacc(:,:,1).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'r--','LineWidth',1);
ploth(5)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_co2_adv_diff,3).*grid.hfacc(:,:,1).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'b','LineWidth',2);
ploth(6)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_surface_fluxes,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'g','LineWidth',2);
plot([-80 80],[0 0],'k--')
set(gca,'XLim',[-80 80],'YLim',[-1e-7 1e-7],'FontSize',12)
xlabel('Latitude','FontSize',12)
ylabel('Weighted zonal average CO2 Flux (mol.m-2.s-1)','FontSize',12)
title('Various CO2 balance residuals','FontSize',12)
plot_names={'Reference','DIC budget','Components ugrad',...
    'Components gradkgrad','Components adv and diff','Surface fluxes'};
legend(plot_names,'Location','SouthWest','Fontsize',12)
orient landscape
eval(['print -dpsc ',obsid,'_zonalwave_co2_flux_components.ps'])
fixpslinestyle([obsid,'_zonalwave_co2_flux_components.ps'])

figure
ploth(1)=plot(grid.latc,squeeze(nanmean(nanmean(dic_co2_flux,3).*grid.hfacc(:,:,1),1)),'k','LineWidth',2);
hold on
ploth(4)=plot(grid.latc,squeeze(nanmean(nanmean(-flux_co2_budget,3).*grid.hfacc(:,:,1),1)),'g--','LineWidth',1);
ploth(2)=plot(grid.latc,squeeze(nanmean(nanmean(-flux_co2_ugradcsat,3).*grid.hfacc(:,:,1),1)),'r','LineWidth',2);
ploth(2)=plot(grid.latc,squeeze(nanmean(nanmean(-flux_co2_gradkgradcsat,3).*grid.hfacc(:,:,1),1)),'r--','LineWidth',1);
ploth(3)=plot(grid.latc,squeeze(nanmean(nanmean(-flux_co2_adv_diff,3).*grid.hfacc(:,:,1),1)),'b','LineWidth',2);
ploth(5)=plot(grid.latc,squeeze(nanmean(nanmean(-flux_surface_fluxes,3).*grid.hfacc(:,:,1),1)),'g','LineWidth',2);
plot([-80 80],[0 0],'k--')
set(gca,'XLim',[-80 80],'YLim',[-1e-7 1e-7],'FontSize',12)
xlabel('Latitude','FontSize',12)
ylabel('Zonally-averaged CO2 Flux (mol.m-2.s-1)','FontSize',12)
title('Various CO2 balance residuals','FontSize',12)
plot_names={'Reference','DIC budget','Components ugrad',...
    'Components gradkgrad','Components adv and diff','Surface fluxes'};
legend(plot_names,'Location','SouthWest','Fontsize',12)
orient landscape
eval(['print -dpsc ',obsid,'_zonalave_co2_flux_components.ps'])
fixpslinestyle([obsid,'_zonalave_co2_flux_components.ps'])

% figure
% ploth(1)=plot(grid.latc,squeeze(nansum(nanmean(dic_co2_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',1);
% hold on
% ploth(2)=plot(grid.latc,squeeze(nansum(nanmean(-flux_co2_components,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',2);
% ploth(3)=plot(grid.latc,squeeze(nansum(nanmean(-ugradtheta.*ddicdth,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'r','LineWidth',2);
% ploth(4)=plot(grid.latc,squeeze(nansum(nanmean(-ugradsalt.*ddicds,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'b','LineWidth',2);
% ploth(5)=plot(grid.latc,squeeze(nansum(nanmean(-ugradalk.*ddicdalk,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g','LineWidth',2);
% ploth(6)=plot(grid.latc,squeeze(nansum(nanmean(-Rcp.*biop,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'r--','LineWidth',1);
% ploth(7)=plot(grid.latc,squeeze(nansum(nanmean(-Rcaco3.*Rcp.*biop,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'b--','LineWidth',1);
% ploth(8)=plot(grid.latc,squeeze(nansum(nanmean(-dic_vrt_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g--','LineWidth',1);
% ploth(9)=plot(grid.latc,squeeze(nansum(nanmean(-(dic_diff_horz+dic_diff_vert),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m','LineWidth',1);
% ploth(10)=plot(grid.latc,squeeze(nansum(nanmean(-ugradcres,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m--','LineWidth',1);
% plot([-80 80],[0 0],'k--')
% set(gca,'XLim',[-80 80],'FontSize',12)
% xlabel('Latitude','FontSize',12)
% ylabel('Weighted zonal average CO2 Flux (mol.m-2.s-1)','FontSize',12)
% title('Components of the "ugrad" CO2 balance','FontSize',12)
% plot_names={'MITgcm CO2 Flux','Total Flux','Csat T Adv','Csat S Adv','Csat A Adv','Soft tissue','Carbonate','DIC Dilution',...
%     'DIC Diffusion','Cres Advection'};
% legend(plot_names,'Location','SouthWest','Fontsize',12)
% orient landscape
% eval(['print -dpsc ',obsid,'_zonal_flux_ugrad_components.ps'])
% fixpslinestyle([obsid,'_zonal_flux_ugrad_components.ps'])
% 
% figure
% ploth(1)=plot(grid.latc,squeeze(nansum(nanmean(dic_co2_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',1);
% hold on
% ploth(2)=plot(grid.latc,squeeze(nansum(nanmean(-flux_co2_ugradcsat,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',2);
% ploth(3)=plot(grid.latc,squeeze(nansum(nanmean(-advtheta.*ddicdth,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'r','LineWidth',2);
% ploth(4)=plot(grid.latc,squeeze(nansum(nanmean(-advsalt.*ddicds,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'b','LineWidth',1);
% ploth(5)=plot(grid.latc,squeeze(nansum(nanmean(-advalk.*ddicdalk,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g','LineWidth',1);
% ploth(6)=plot(grid.latc,squeeze(nansum(nanmean(-difftheta.*ddicdth,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'r--','LineWidth',2);
% ploth(7)=plot(grid.latc,squeeze(nansum(nanmean(-diffsalt.*ddicds,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'b--','LineWidth',1);
% ploth(8)=plot(grid.latc,squeeze(nansum(nanmean(-diffalk.*ddicdalk,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g--','LineWidth',1);
% ploth(9)=plot(grid.latc,squeeze(nansum(nanmean(-Rcp.*biop,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m-.','LineWidth',1);
% ploth(10)=plot(grid.latc,squeeze(nansum(nanmean(Rcp.*(po4_adv_horz+po4_adv_vert),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'c','LineWidth',1);
% ploth(11)=plot(grid.latc,squeeze(nansum(nanmean(Rcp.*(po4_diff_horz+po4_diff_vert),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'c--','LineWidth',1);
% ploth(12)=plot(grid.latc,squeeze(nansum(nanmean(-Rcaco3.*Rcp.*biop+(Rnp.*biop),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m:','LineWidth',1);
% ploth(13)=plot(grid.latc,squeeze(nansum(nanmean(-dic_vrt_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g:','LineWidth',1);
% ploth(14)=plot(grid.latc,squeeze(nansum(nanmean(-advcres,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m','LineWidth',1);
% ploth(15)=plot(grid.latc,squeeze(nansum(nanmean(-diffcres,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m--','LineWidth',1);
% ploth(16)=plot(grid.latc,squeeze(nansum(nanmean(-diffcres-advcres-Rcp.*biop-Rcaco3.*Rcp.*biop+(Rnp.*biop),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'m','LineWidth',2);
% %fw_components=(advsalt+diffsalt).*ddicds+(advalk+diffalk).*ddicdalk+dic_vrt_flux;
% % Should subtract biological effect on alkalinity here
% fw_components=(advsalt+diffsalt).*ddicds+(advalk+diffalk-(Rnp.*biop)+(Rcaco3.*Rcp.*biop)).*ddicdalk+dic_vrt_flux;
% %fw_components=(advsalt+diffsalt).*ddicds+(advalk+diffalk+alk_carb_flux+alk_bio_flux).*ddicdalk+dic_vrt_flux;
% ploth(17)=plot(grid.latc,squeeze(nansum(nanmean(-fw_components,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'g','LineWidth',2);
% plot([-80 80],[0 0],'k--')
% set(gca,'XLim',[-80 80],'YLim',[-1e-7 1e-7],'FontSize',12)
% xlabel('Latitude','FontSize',12)
% ylabel('Weighted zonal average CO2 Flux (mol.m-2.s-1)','FontSize',12)
% title('Components of the "advection and diffusion" CO2 balance','FontSize',12)
% plot_names={'MITgcm CO2 Flux','Total Flux','Csat T Adv','Csat S Adv','Csat A Adv',...
%     'Csat T Diff','Csat S Diff','Csat A Diff','Soft tissue','P Adv','P Diff','Carbonate','DIC Dilution',...
%     'Cres Advection','Cres Diffusion','Cres+Bio+Carb','Salt+Alk+FW'};
% legend(plot_names,'Location','SouthWest','Fontsize',12)
% orient landscape
% eval(['print -dpsc ',obsid,'_zonal_flux_adv_diff_components.ps'])
% fixpslinestyle([obsid,'_zonal_flux_adv_diff_components.ps'])
%
figure
ploth(1)=plot(grid.latc,squeeze(nansum(nanmean(dic_co2_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),'k','LineWidth',1);
hold on
ploth(2)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-flux_surface_fluxes,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'k','LineWidth',2);
ploth(3)=plot(grid.latc,smooth(squeeze(nansum(nanmean(theta_surf_flux.*ddicdth,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'r','LineWidth',2);
ploth(4)=plot(grid.latc,smooth(squeeze(nansum(nanmean(salt_surf_flux.*ddicds,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'b','LineWidth',2);
ploth(5)=plot(grid.latc,smooth(squeeze(nansum(nanmean(alk_vrt_flux.*ddicdalk,3).*grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)./nansum(grid.dxc(:,:,1),1)),11,'rlowess'),'g:','LineWidth',1);
ploth(6)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-(Rcp.*biop),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'m-.','LineWidth',1);
ploth(7)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-Rcaco3.*Rcp.*biop+(Rnp.*biop),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'m:','LineWidth',1);
ploth(8)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-dic_vrt_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'g--','LineWidth',1);
ploth(9)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-cres_surface_fluxes,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'m--','LineWidth',1);
ploth(10)=plot(grid.latc,smooth(squeeze(nansum(nanmean(-cres_surface_fluxes-(Rcp.*biop)-Rcaco3.*Rcp.*biop+(Rnp.*biop),3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'m','LineWidth',2);
ploth(11)=plot(grid.latc,smooth(squeeze(nansum(nanmean(salt_surf_flux.*ddicds+alk_vrt_flux.*ddicdalk-dic_vrt_flux,3).*grid.dxc(:,:,1),1)./nansum(grid.dxc(:,:,1).*grid.hfacc(:,:,1),1)),11,'rlowess'),'g','LineWidth',2);
plot([-80 80],[0 0],'k--')
set(gca,'XLim',[-80 80],'YLim',[-2e-7 2e-7],'FontSize',12)
xlabel('Latitude','FontSize',12)
ylabel('Weighted zonal average CO2 Flux (mol.m-2.s-1)','FontSize',12)
title('Components of the surface flux CO2 balance','FontSize',12)
plot_names={'MITgcm CO2 Flux','Total Flux','Csat T flux','Csat S flux','Csat A flux','Soft tissue','Carbonate','DIC Dilution',...
    'Cres Advdiff','Cres+Bio+Carb','Salt+Alk+FW'};
legend(plot_names,'Location','SouthWest','Fontsize',12)
orient landscape
eval(['print -dpsc ',obsid,'_zonal_surface_flux_components.ps'])
fixpslinestyle([obsid,'_zonal_surface_flux_components.ps'])

%% Components of the "flux_co2_ugradcsat" breakdown

% figure
% h(1)=subplot(331);
% [~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(-ugradtheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-ugradtheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Csat T';'Transport [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(2)=subplot(332);
% [~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(-ugradsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-ugradsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Csat S';'Transport [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(3)=subplot(333);
% [~,a(3)]=contourf(grid.lonc,grid.latc,nanmean(-ugradalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-ugradalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Csat A';'Transport [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(4)=subplot(334);
% [~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(-ugradcres,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-ugradcres,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Transport [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(5)=subplot(335);
% [~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(-((Rcp.*biop)+(Rcaco3.*Rcp.*biop)),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-((Rcp.*biop)+(Rcaco3.*Rcp.*biop)),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to PP and';'Carbonate [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(6)=subplot(336);
% [~,a(6)]=contourf(grid.lonc,grid.latc,nanmean(-(dic_diff_vert+dic_diff_horz),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-(dic_diff_vert+dic_diff_horz),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to DIC';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(7)=subplot(337);
% [~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'Surface Dilution of DIC';'[mol.m-2.s-1]'},'FontSize',12);xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(8)=subplot(338);
% [~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(-ugradsalt.*ddicds-ugradalk.*ddicdalk-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-ugradsalt.*ddicds-ugradalk.*ddicdalk-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Fluxes due to sum of';'FW Components [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(9)=subplot(339);
% [~,a(9)]=contourf(grid.lonc,grid.latc,nanmean(-flux_co2_ugradcsat,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-flux_co2_ugradcsat,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Fluxes due to';'Components [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% set([a;b],'LineStyle','none')
% set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
% orient landscape
% eval(['print -dpsc ',obsid,'_flux_co2_ugradcsat.ps'])

%% The breakdown of cresidual into horizontal and vertical components
% figure
% h(1)=subplot(321);
% [~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(-advcres,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-advcres,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Total Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(2)=subplot(322);
% [~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(-diffcres,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-diffcres,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Total Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(3)=subplot(323);
% [~,a(3)]=contourf(grid.lonc,grid.latc,nanmean(-cres_adv_horz,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(cres_adv_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Horizontal Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(4)=subplot(324);
% [~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(-cres_diff_horz,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-cres_diff_horz,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Horizontal Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(5)=subplot(325);
% [~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(-cres_adv_vert,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-cres_adv_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Vertical Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% h(6)=subplot(326);
% [~,a(6)]=contourf(grid.lonc,grid.latc,nanmean(-cres_diff_vert,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
% tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
% hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-cres_diff_vert,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
% b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
% title({'CO2 Flux due to Cres';'Vertical Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
% set([a;b],'LineStyle','none')
% set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
% orient landscape
% eval(['print -dpsc ',obsid,'_cres_fluxes.ps'])

%% The flux_surface_fluxes breakdown
figure
h(1)=subplot(331);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(theta_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Csat';'Heat flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(332);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Csat';'Salt flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(333);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean(alk_vrt_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_vrt_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Csat';'Alk flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(334);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to DIC';'Dilution [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(335);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(-cres_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-cres_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Cres';'Advection/Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(336);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean(-Rcp.*biop-Rcaco3.*Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-Rcp.*biop-Rcaco3.*Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Net Biological';'Activity [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(337);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(-cres_surface_fluxes-Rcp.*biop-Rcaco3.*Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-cres_surface_fluxes-Rcp.*biop-Rcaco3.*Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to';'Cres+Bio [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(338);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds+alk_vrt_flux.*ddicdalk-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds+alk_vrt_flux.*ddicdalk-dic_vrt_flux,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Fluxes due to sum of';'FW Components [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(339);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean(-flux_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-flux_surface_fluxes,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Fluxes due to sum of';'All Components [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_flux_surface_fluxes.ps'])

%% Plot the breakdown between fluxes caused by surface forcing and relaxation
% if strcmp(obsid(1:6),'mitgcm')
%     trelax_surf_flux=(squeeze(surfdiag.TRELAX).*repmat(grid.hfacc(:,:,1),[1,1,size(surfdiag.TFLUX,4)]))./(rhoconst*cp); % K m s-1
%     tforc_surf_flux=theta_surf_flux-trelax_surf_flux;
%     srelax_surf_flux=squeeze(surfdiag.SRELAX).*repmat(grid.hfacc(:,:,1),[1,1,size(surfdiag.SRELAX,4)])./rhoconst;
%     sforc_surf_flux=salt_surf_flux-srelax_surf_flux;
%     
%     surf_salt = salt(:,:,1); smean=nansum(surf_salt(:).*grid.rac(:))./nansum(grid.rac(:));
%     surf_alk  = alk(:,:,1);  amean=nansum(surf_alk (:).*grid.rac(:))./nansum(grid.rac(:));
%     alk_vrtrelax_flux=srelax_surf_flux.*(amean/smean);
%     alk_vrtforc_flux=sforc_surf_flux.*(amean/smean);
%     
%     figure
%     h(1)=subplot(331);
%     [~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(theta_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(theta_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Heat flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(2)=subplot(332);
%     [~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(salt_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Salt flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(3)=subplot(333);
%     [~,a(3)]=contourf(grid.lonc,grid.latc,nanmean(alk_vrt_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_vrt_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
% %    hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean((alk_vrt_flux-Rnp.*biop+Rcaco3.*Rcp.*biop).*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Alk flux [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(4)=subplot(334);
%     [~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(tforc_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(tforc_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Heat Forcing [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(5)=subplot(335);
%     [~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(sforc_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(sforc_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Salt Forcing [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(6)=subplot(336);
%     [~,a(6)]=contourf(grid.lonc,grid.latc,nanmean(alk_vrtforc_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_vrtforc_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Alk (FW) Forcing [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(7)=subplot(337);
%     [~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(trelax_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(trelax_surf_flux.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Heat Relaxation [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(8)=subplot(338);
%     [~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(srelax_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(srelax_surf_flux.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Salt Relaxation [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     h(9)=subplot(339);
%     [~,a(9)]=contourf(grid.lonc,grid.latc,nanmean(alk_vrtrelax_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(alk_vrtrelax_flux.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end;
%     b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop
%     title({'CO2 Flux due to Csat';'Alk (FW) Relaxation [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
%     set([a;b],'LineStyle','none')
%     set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12);
%     orient landscape
%     eval(['print -dpsc ',obsid,'_flux_forcrelax_fluxes.ps'])
% end

%% flux_co2_adv_diff
figure
h(1)=subplot(341);
[~,a(1)]=contourf(grid.lonc,grid.latc,nanmean(-advtheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(1),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-advtheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(1)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to T';'Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(2)=subplot(342);
[~,a(2)]=contourf(grid.lonc,grid.latc,nanmean(-advsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(2),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-advsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(2)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to S';'Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(3)=subplot(343);
[~,a(3)]=contourf(grid.lonc,grid.latc,nanmean(-advalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(3),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-advalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(3)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to A';'Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(4)=subplot(344);
[~,a(4)]=contourf(grid.lonc,grid.latc,nanmean(-advcres,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(4),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-advcres,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(4)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Cres';'Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(5)=subplot(345);
[~,a(5)]=contourf(grid.lonc,grid.latc,nanmean(-difftheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(5),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-difftheta.*ddicdth,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(5)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to T';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(6)=subplot(346);
[~,a(6)]=contourf(grid.lonc,grid.latc,nanmean(-diffsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(6),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-diffsalt.*ddicds,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(6)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to S';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(7)=subplot(347);
[~,a(7)]=contourf(grid.lonc,grid.latc,nanmean(-diffalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(7),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-diffalk.*ddicdalk,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(7)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to A';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(8)=subplot(348);
[~,a(8)]=contourf(grid.lonc,grid.latc,nanmean(-diffcres,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(8),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-diffcres,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(8)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Cres';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(9)=subplot(349);
[~,a(9)]=contourf(grid.lonc,grid.latc,nanmean(Rcp.*(po4_adv_horz+po4_adv_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(9),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(Rcp.*(po4_adv_horz+po4_adv_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(9)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Phosphate';'Advection [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(10)=subplot(3,4,10);
[~,a(10)]=contourf(grid.lonc,grid.latc,nanmean(Rcp.*(po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(10),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(Rcp.*(po4_diff_horz+po4_diff_vert),3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(10)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Flux due to Phosphate';'Diffusion [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(11)=subplot(3,4,11);
[~,a(11)]=contourf(grid.lonc,grid.latc,nanmean(-diffcres-advcres-Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(11),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-diffcres-advcres-Rcp.*biop,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(11)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Fluxes due to';'Cres+Bio [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
h(12)=subplot(3,4,12);
fw_components=(advsalt+diffsalt).*ddicds+(advalk+diffalk).*ddicdalk+dic_vrt_flux;
[~,a(12)]=contourf(grid.lonc,grid.latc,nanmean(-fw_components,3)'.*grid.hfacc(:,:,1)',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
tmp=get(a(12),'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
hold on;  if isempty(strfind(obsid,'ecco2darwin')); contour(grid.lonc,grid.latc,nanmean(-fw_components,3)'.*grid.hfacc(:,:,1)',[0 0],'k'); end; 
if strcmpi(plot_landmask,'true'); b(12)=jpcolor(grid.lonc,grid.latc,landmask.*-2e-7); axestop; end
title({'CO2 Fluxes due to sum of';'FW Components [mol.m-2.s-1]'});xlabel('Longitude','FontSize',12);ylabel('Latitude','FontSize',12)
set([a;b],'LineStyle','none')
set(h,'Xlim',[min(grid.long) max(grid.lonc)],'XTickLabelMode','manual','XTickMode','manual','XTick',[min(grid.lonc):120:max(grid.long)],'XTickLabel',[floor(min(grid.long)):120:ceil(max(grid.lonc))],'YLim',[min(grid.latc) max(grid.latc)],'YTick',[-80:40:80],'FontSize',12); 
orient landscape
eval(['print -dpsc ',obsid,'_flux_co2_adv_diff.ps'])
clear a b h 

%% Barchart of areal integrations
if grid.nx==128 && grid.ny==64
    cflux.eqst.mask=grid.hfacc(:,:,1);
    cflux.eqst.mask(:,[1:16,49:64])=NaN;
    cflux.eqst.mask(10:128,1:19)=NaN; %extra bit in Southern Ocean
    cflux.eqst.mask(nanmean(dic_co2_flux,3)>0)=NaN;
    %eqst_lat=nansum(eqst_region(:).*grid.yc(:).*grid.rac(:))./nansum(eqst_region(:).*grid.rac(:));
    cflux.eqst.mask(111:120,36:40)=1; %extra bit in North Atlantic
    cflux.eqst.mask(75:77,47:49)=NaN; %extra bit in North Atlantic
    cflux.eqst.lat=0;
        
    cflux.aaiw.mask=grid.hfacc(:,:,1);
    cflux.aaiw.mask(:,[1:12,32:64])=NaN;
    cflux.aaiw.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.aaiw.mask(nanmean(dic_co2_flux,3)<0)=NaN;
    cflux.aaiw.mask(105:128,15:20)=1; % Add South Atlantic back in
    %aaiw_lat=nansum(aaiw_region(:).*grid.yc(:).*grid.rac(:))./nansum(aaiw_region(:).*grid.rac(:));
    cflux.aaiw.lat=-40;
    
    cflux.aabw.mask=grid.hfacc(:,:,1);
    cflux.aabw.mask(:,32:64)=NaN;
    cflux.aabw.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.aabw.mask(~isnan(cflux.aaiw.mask))=NaN;
    %aabw_lat=nansum(aabw_region(:).*grid.yc(:).*grid.rac(:))./nansum(aabw_region(:).*grid.rac(:));
    cflux.aabw.lat=-65;
    
    cflux.npna.mask=grid.hfacc(:,:,1);
    cflux.npna.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.npna.mask(~isnan(cflux.aaiw.mask))=NaN;
    cflux.npna.mask(~isnan(cflux.aabw.mask))=NaN;
    %npna_lat=nansum(npna_region(:).*grid.yc(:).*grid.rac(:))./nansum(npna_region(:).*grid.rac(:));
    cflux.npna.lat=45;

    %global
    cflux.glob.mask=grid.hfacc(:,:,1);
    cflux.glob.lat=0;
    cflux.glob.lon=180;    

    % North West Atlantic
    cflux.nwa.mask=cflux.npna.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.nwa.mask(119:128,:)=NaN;
    cflux.nwa.lat=45;
    cflux.nwa.lon=300;
    % North East Atlantic
    cflux.nea.mask=cflux.npna.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.nea.mask(64:118,:)=NaN;
    cflux.nea.lat=45;
    cflux.nea.lon=345;
    % Equatorial/subtropical Atlantic
    cflux.ea.mask=cflux.eqst.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.ea.lat=0;
    cflux.ea.lon=330;
    % NW Pacific
    cflux.nwp.mask=cflux.npna.mask.*grid.pacific_hfacc(:,:,1);
    cflux.nwp.mask(65:128,:)=NaN;
    cflux.nwp.lat=45;
    cflux.nwp.lon=170;
    % NE Pacific
    cflux.nep.mask=cflux.npna.mask.*grid.pacific_hfacc(:,:,1);
    cflux.nep.mask(1:64,:)=NaN;
    cflux.nep.lat=45;
    cflux.nep.lon=240;
    % W Equatorial Pacific
    cflux.wep.mask=cflux.eqst.mask.*grid.pacific_hfacc(:,:,1);
    cflux.wep.mask(65:128,:)=NaN;
%    mitgcm.wep.mask(42:46,26:30)=NaN;
    cflux.wep.lat=0;
    cflux.wep.lon=170;
    % E Equatorial Pacific
    cflux.eep.mask=cflux.eqst.mask.*grid.pacific_hfacc(:,:,1);
    cflux.eep.mask(1:64,:)=NaN;
    cflux.eep.lat=0;
    cflux.eep.lon=240;
    % Equatorial Indian Ocean
    cflux.ei.mask=cflux.eqst.mask.*grid.indic_hfacc(:,:,1);
    cflux.ei.lat=0;
    cflux.ei.lon=75;
    % SO Indian Sector 
    cflux.soi.mask=cflux.aaiw.mask.*grid.indic_hfacc(:,:,1);
    cflux.soi.lat=-40;
    cflux.soi.lon=75;
    % SO Pacific Sector
    cflux.sop.mask=cflux.aaiw.mask.*grid.pacific_hfacc(:,:,1);
%    mitgcm.sop.mask(42:46,1:21)=NaN;
    cflux.sop.lat=-40;
    cflux.sop.lon=190;
    % SO Atlantic Sector
    cflux.soa.mask=cflux.aaiw.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.soa.lat=-40;
    cflux.soa.lon=330;
    % East Antarctica (Indian)
    cflux.eai.mask=cflux.aabw.mask.*grid.indic_hfacc(:,:,1);
    cflux.eai.lat=-65;
    cflux.eai.lon=75;
    % Ross Sea (Pacific)
    cflux.rsp.mask=cflux.aabw.mask.*grid.pacific_hfacc(:,:,1);
    cflux.rsp.lat=-65;
    cflux.rsp.lon=190;
    % Weddell Sea (Atlantic)
    cflux.wsa.mask=cflux.aabw.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.wsa.lat=-65;
    cflux.wsa.lon=330;
elseif grid.nx==360 && grid.ny>=160    
%%%%%%%%%%%
    cflux.eqst.mask=grid.hfacc(:,:,1);
    cflux.eqst.mask(:,[1:35,127:160])=NaN;
    cflux.eqst.mask(nanmean(dic_co2_flux,3)>0)=NaN;
    %eqst_lat=nansum(eqst_region(:).*grid.yc(:).*grid.rac(:))./nansum(eqst_region(:).*grid.rac(:));
    cflux.eqst.lat=0;
        
    cflux.aaiw.mask=grid.hfacc(:,:,1);
    cflux.aaiw.mask(:,[1:23,80:158])=NaN;
    cflux.aaiw.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.aaiw.mask(nanmean(dic_co2_flux,3)<0)=NaN;
    cflux.aaiw.lat=-40;
    
    cflux.aabw.mask=grid.hfacc(:,:,1);
    cflux.aabw.mask(:,80:160)=NaN;
    cflux.aabw.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.aabw.mask(~isnan(cflux.aaiw.mask))=NaN;
    cflux.aabw.mask(260:275,30:40)=NaN;
    cflux.aabw.lat=-65;
    
    cflux.npna.mask=grid.hfacc(:,:,1);
    cflux.npna.mask(~isnan(cflux.eqst.mask))=NaN;
    cflux.npna.mask(~isnan(cflux.aaiw.mask))=NaN;
    cflux.npna.mask(~isnan(cflux.aabw.mask))=NaN;
    cflux.npna.mask(260:275,30:40)=NaN;
    cflux.npna.lat=45;
    
    %global
    cflux.glob.mask=grid.hfacc(:,:,1);
    cflux.glob.lat=0;
    cflux.glob.lon=180;
    
    % North West Atlantic
    cflux.nwa.mask=cflux.npna.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.nwa.mask(334:360,:)=NaN;
    cflux.nwa.mask(1:70,:)=NaN;  
    cflux.nwa.lat=45;
    cflux.nwa.lon=300;
    % North East Atlantic
    cflux.nea.mask=cflux.npna.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.nea.mask(180:331,:)=NaN;
    cflux.nea.lat=45;
    cflux.nea.lon=345;
    % Equatorial/subtropical Atlantic
    cflux.ea.mask=cflux.eqst.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.ea.lat=0;
    cflux.ea.lon=330;
    % NW Pacific
    cflux.nwp.mask=cflux.npna.mask.*grid.pacific_hfacc(:,:,1);
    cflux.nwp.mask(181:360,:)=NaN;
    cflux.nwp.lat=45;
    cflux.nwp.lon=170;
    % NE Pacific
    cflux.nep.mask=cflux.npna.mask.*grid.pacific_hfacc(:,:,1);
    cflux.nep.mask(1:180,:)=NaN;
    cflux.nep.lat=45;
    cflux.nep.lon=240;
    % W Equatorial Pacific
    cflux.wep.mask=cflux.eqst.mask.*grid.pacific_hfacc(:,:,1);
    cflux.wep.mask(181:360,:)=NaN;
    cflux.wep.lat=0;
    cflux.wep.lon=170;
    % E Equatorial Pacific
    cflux.eep.mask=cflux.eqst.mask.*grid.pacific_hfacc(:,:,1);
    cflux.eep.mask(1:180,:)=NaN;
    cflux.eep.lat=0;
    cflux.eep.lon=240;
    % Equatorial Indian Ocean
    cflux.ei.mask=cflux.eqst.mask.*grid.indic_hfacc(:,:,1);
    cflux.ei.lat=0;
    cflux.ei.lon=75;
    % SO Indian Sector 
    cflux.soi.mask=cflux.aaiw.mask.*grid.indic_hfacc(:,:,1);
    cflux.soi.lat=-40;
    cflux.soi.lon=75;
    % SO Pacific Sector
    cflux.sop.mask=cflux.aaiw.mask.*grid.pacific_hfacc(:,:,1);
    cflux.sop.lat=-40;
    cflux.sop.lon=190;
    % SO Atlantic Sector
    cflux.soa.mask=cflux.aaiw.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.soa.lat=-40;
    cflux.soa.lon=330;
    % East Antarctica (Indian)
    cflux.eai.mask=cflux.aabw.mask.*grid.indic_hfacc(:,:,1);
    cflux.eai.lat=-65;
    cflux.eai.lon=75;
    % Ross Sea (Pacific)
    cflux.rsp.mask=cflux.aabw.mask.*grid.pacific_hfacc(:,:,1);
    cflux.rsp.lat=-65;
    cflux.rsp.lon=190;
    % Weddell Sea (Atlantic)
    cflux.wsa.mask=cflux.aabw.mask.*grid.atlantic_hfacc(:,:,1);
    cflux.wsa.lat=-65;
    cflux.wsa.lon=330; 
end
%%
    % Calculate weighted average fluxes over the region
    region={'glob','npna','aabw','aaiw','eqst','nwa','nea','ea','nwp','nep','wep','eep','ei','soi','sop','soa','eai','rsp','wsa'};
    component={...
        'dic_co2_flux',...
        'theta_surf_flux.*ddicdth',...
        'salt_surf_flux.*ddicds',...
        'alk_vrt_flux.*ddicdalk',...
        '-dic_vrt_flux',...
        '-cres_surface_fluxes',...
        '-Rcp.*biop',...
        '-Rcaco3.*Rcp.*biop',...
        '-cres_surface_fluxes-Rcp.*biop',...
        'salt_surf_flux.*ddicds+alk_vrt_flux.*ddicdalk-dic_vrt_flux',...
        '-advtheta.*ddicdth',...
        '-advsalt.*ddicds',...
        '-advalk.*ddicdalk',...
        '-advcres',...
        '-difftheta.*ddicdth',...
        '-diffsalt.*ddicds',...
        '-diffalk.*ddicdalk',...
        '-diffcres',...
        'Rcp.*(po4_adv_horz+po4_adv_vert)',...
        'Rcp.*(po4_diff_horz+po4_diff_vert)',...
        '-diffcres-advcres-Rcp.*biop',...
        '-fw_components',...
        'dic_adv_horz+dic_adv_vert',...
        'dic_diff_horz+dic_diff_vert',...
        'dic_bio_flux+dic_carb_flux',...
        'dic_vrt_flux',...
        '-flux_surface_fluxes',...
        };
    
    names={...
        'dic_co2_flux',...
        'theta_surf_flux',...
        'salt_surf_flux',...
        'alk_vrt_flux',...
        'dic_vrt_flux',...
        'cres_co2_flux',...
        'bio_co2_flux',...
        'carb_co2_flux',...
        'net_bio_flux',...
        'net_fw_flux',...
        'theta_adv_flux',...
        'salt_adv_flux',...
        'alk_adv_flux',...
        'cres_adv_flux',...
        'theta_diff_flux',...
        'salt_diff_flux',...
        'alk_diff_flux',...
        'cres_diff_flux',...
        'po4_adv_flux',...
        'po4_diff_flux',...
        'net_bio_advdiff',...
        'net_fw_advdiff',...
        'dic_advection,',...
        'dic_diffusion',...
        'dic_bio_activity',...
        'dic_vrt_flux',...
        'flux_surface_fluxes',...
        };
    
    for i=1:length(region)
        cflux.(region{i}).area=nansum(grid.rac(:).*cflux.(region{i}).mask(:));
        for j=1:length(component)
            eval(['tmp=nanmean(',component{j},',3);'])
            cflux.(region{i}).components(j)=nansum(tmp(:).*grid.rac(:).*cflux.(region{i}).mask(:)).*12e-15*(60*60*24*360); % PgC/yr ./nansum(grid.rac(:).*mitgcm.(region{i}).mask(:)); % mol m-2 s-1
        end
    end
    
    % create cmap consistent with bluewhitered colormap for bars
    bottom = [0 0 0.5];
    botmiddle = [0 0.5 1];
    middle = [1 1 1];
    topmiddle = [1 0 0];
    top = [0.5 0 0];
    newmap = zeros(12, 3);
    cmap = [bottom; botmiddle; middle; topmiddle; top];
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(linspace(0, 1, 5), cmap(:,i), linspace(0, 1, 12))', 0), 1);
    end
       
    % 1.) Net totals
    figure
    subplot(221)
    tot_vec=[1,27,2,9,10];
    bar([cflux.npna.lat,cflux.eqst.lat,cflux.aaiw.lat,cflux.aabw.lat],...
    [cflux.npna.components(tot_vec);cflux.eqst.components(tot_vec);cflux.aaiw.components(tot_vec);cflux.aabw.components(tot_vec)])
    colormap(bluewhitered(3));legend(names(tot_vec),'Interpreter','none')
    set(gca,'XLim',[-80 80],'YLim',[-2.5 2.5],'XTickLabelMode','manual','XTickMode','manual','XTick',[-80:20:80],'XTickLabel',[-80:20:80])
    xlabel('Latitude','FontSize',12);ylabel('Integrated CO2 fluxes [PgC/yr]','FontSize',12)
    
    % 2.) Surface Fluxes
    subplot(222)
    surf_vec=[2:8];
    bar([cflux.npna.lat,cflux.eqst.lat,cflux.aaiw.lat,cflux.aabw.lat],...
    [cflux.npna.components(surf_vec);cflux.eqst.components(surf_vec);cflux.aaiw.components(surf_vec);cflux.aabw.components(surf_vec)])
    colormap(bluewhitered(7));legend(names(surf_vec),'Interpreter','none')
    set(gca,'XLim',[-80 80],'YLim',[-2.5 2.5],'XTickLabelMode','manual','XTickMode','manual','XTick',[-80:20:80],'XTickLabel',[-80:20:80])
    xlabel('Latitude','FontSize',12);ylabel('Integrated CO2 fluxes [PgC/yr]','FontSize',12)

    % 3.) Net Totals again
    subplot(223)
    tot_vec2=[1,27,2,21,22];
    bar([cflux.npna.lat,cflux.eqst.lat,cflux.aaiw.lat,cflux.aabw.lat],...
    [cflux.npna.components(tot_vec2);cflux.eqst.components(tot_vec2);cflux.aaiw.components(tot_vec2);cflux.aabw.components(tot_vec2)])
    colormap(bluewhitered(3));legend(names(tot_vec2),'Interpreter','none')
    set(gca,'XLim',[-80 80],'YLim',[-2.5 2.5],'XTickLabelMode','manual','XTickMode','manual','XTick',[-80:20:80],'XTickLabel',[-80:20:80])
    xlabel('Latitude','FontSize',12);ylabel('Integrated CO2 fluxes [PgC/yr]','FontSize',12)
    
    % 4.) Advective and diffusive fluxes
    subplot(224)
    ad_vec=[11,15,12,16,13,17,5,14,18,19,20,8];
    bar([cflux.npna.lat,cflux.eqst.lat,cflux.aaiw.lat,cflux.aabw.lat],...
    [cflux.npna.components(ad_vec);cflux.eqst.components(ad_vec);cflux.aaiw.components(ad_vec);cflux.aabw.components(ad_vec)])
    colormap(newmap);legend(names(ad_vec),'Interpreter','none')
    set(gca,'XLim',[-80 80],'YLim',[-2.5 2.5],'XTickLabelMode','manual','XTickMode','manual','XTick',[-80:20:80],'XTickLabel',[-80:20:80])
    xlabel('Latitude','FontSize',12);ylabel('Integrated CO2 fluxes [PgC/yr]','FontSize',12)
    orient landscape
    eval(['print -dpsc ',obsid,'_regional_bar_charts.ps'])
   
    
    bar_vec=[27,2,7,6,10];
    % Write some data out for GMT
    dlmwrite('co2flux.dat',[]) % Open and clobber a blank file
%    adj=[-7.5,-2.5,2.5,7.5];
    adj=[-10,-5,0,5,10];
    for i=6:length(region)
        for j=1:length(bar_vec)
           dlmwrite('co2flux.dat',[cflux.(region{i}).lon+adj(j),cflux.(region{i}).lat,cflux.(region{i}).components(bar_vec(j)),(100/(length(bar_vec)-1)).*(j-1)],'-append')
        end
    end
    
    % 4.) DIC Budget Totals
    figure
    tot_vec=[1,23,24,25,26];
    adj=[-5,-2.5,0,2.5,5];
    xloc=[repmat(cflux.npna.lat,[1,length(tot_vec)])+adj;repmat(cflux.eqst.lat,[1,length(tot_vec)])+adj;...
        repmat(cflux.aaiw.lat,[1,length(tot_vec)])+adj;repmat(cflux.aabw.lat,[1,length(tot_vec)])+adj];
    bar(xloc,[cflux.npna.components(tot_vec);cflux.eqst.components(tot_vec);cflux.aaiw.components(tot_vec);cflux.aabw.components(tot_vec)])
    colormap(bluewhitered(5));caxis([1 5]);legend(names(tot_vec),'Interpreter','none')
    set(gca,'XLim',[-80 80],'YLim',[-2.5 2.5],'XTickLabelMode','manual','XTickMode','manual','XTick',[-80:20:80],'XTickLabel',[-80:20:80])
    xlabel('Latitude','FontSize',12);ylabel('Integrated CO2 fluxes [PgC/yr]','FontSize',12)
    orient landscape
    eval(['print -dpsc ',obsid,'_dicbudget_bar_charts.ps'])
% else
%     %Use accumarray to bin into takahashi grid
%     
%     takadx=unique(diff(unique(clim_pco2_lonc)));
%     takady=unique(diff(unique(clim_pco2_latc)));
%     
%     [colID, rowID] = meshgrid( ceil(grid.lonc/takadx), ceil(grid.latc/takady)-min(ceil(grid.latc/takady))+1 );
%     
%     ind=transpose(colID + (rowID-1) .* length(unique(clim_pco2_lonc))) ;
%     flux_ann=nanmean(flux_surface_fluxes,3);
%     flux_binned  = accumarray(ind(:),flux_ann(:).*grid.rac(:)) ./ accumarray(ind(:),grid.rac(:)) ;
%     flux_binned=reshape(flux_binned,[length(unique(clim_pco2_lonc)), length(unique(clim_pco2_latc))]);
%     
%     figure
%     [~,z]=contourf(unique(clim_pco2_lonc),unique(clim_pco2_latc),-flux_binned',[-2e-6,-2e-7:2e-8:2e-7,2e-6]);
%     tmp=get(z,'LevelList');caxis([tmp(2) tmp(end-1)]);colormap(bluewhitered(length(tmp)-3));colorbar('YTick',tmp(2:2:end-1));
%     clear z
%  end
clear a b h
