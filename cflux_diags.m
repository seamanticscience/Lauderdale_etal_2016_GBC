% Compute CO2 flux component diagnostics
% Declare global variables
global cp KDOPRemin seconds_per_year rhoconst kave grid mld_array ntsteps
global perl_2_perkg perm3_2_perl molperm3_2_umolperkg umolperkgmpers_2_molperm2peryr
global Rcn Rcp Rpo Rno Rnp Rcaco3 Rsip
global obsid vertmix ent_method eddy_method kave_method modeldir

% Choose a flavour
obsid='mitgcm_fluxes' % use model flux fields
%obsid='mitgcm_monthly' % use model flux fields

%obsid='mitgcm_monthly_kpp'
%obsid='mitgcm_annual_kpp_mld'

% How to calculate entrainment flux
vertmix=''; % ie use model fluxes
%vertmix='diffusion' % ie calculate diffusive fluxes

% Do you want Eulerian-mean or residual velocities (ie include GM)?
eddy_method='residual'
%eddy_method='eulerian';

% kave method
kave_method='arb_depth'
kave=55 % arbitrary depth in metres
%kave_method='maxmld'
%kave_method='minmld'
%kave_method='meanmld'
%kave_method='mld'

% Set advection scheme to simple center-second orger.
advscheme   =2;
ptradvscheme=77;

%some conversions
cp = 3985.0; % Heat capacity J kg-1 K-1
KDOPRemin=1/(6*30*86400); % six months in s-1
seconds_per_year=360.*86400;
rhoconst=1024.5;
perl_2_perkg=1000/rhoconst; % e.g. mol l-1 * 1000 l m-3 * 1/rho m3 kg-1 -> mol kg-1
perm3_2_perl = 1.0/1000; % e.g. mol m-3 * 1/1000 m3 l-1 -> mol l-1
molperm3_2_umolperkg=1e6/rhoconst;
umolperkgmpers_2_molperm2peryr = rhoconst*seconds_per_year*1e-6 ;

%% set Redfield ratios for organic matter
    Rcn=117/16;
    Rcp=117;
    Rpo=-1/170;
    Rno=-16/170;
    Rnp=16/1;
    Rcaco3=7e-2; % Inorganic/organic carbon rain ratio
    Rsip  =15;
    
    %% Use mitgcm data
    usemld='model'
    %modeldir='';
    modeldir='model_kpp'
       
    start_iter=mit_getparm([modeldir,'/data'],'nIter0');

    % Load model domain
    grid=mit_loadgrid(modeldir);
    grid.dzc=[grid.zc(1);diff(grid.zc)];
    
    %get ocean basin masks
    grid=mit_oceanmasks(grid);
    
    grid.cmask(isnan(grid.cmask))=0;
    grid.umask(isnan(grid.umask))=0;
    grid.vmask(isnan(grid.vmask))=0;
    
    tave=rdmnc([modeldir,'/tave.',num2str(start_iter,'%010d'),'.glob.nc']);
    surfdiag=rdmnc([modeldir,'/surfDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    ocediag=rdmnc([modeldir,'/oceDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    flxdiag=rdmnc([modeldir,'/flxDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    dicflx=rdmnc([modeldir,'/dic_flxDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    alkflx=rdmnc([modeldir,'/alk_flxDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    po4flx=rdmnc([modeldir,'/po4_flxDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    dopflx=rdmnc([modeldir,'/dop_flxDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    ptracers=rdmnc([modeldir,'/ptr_tave.',num2str(start_iter,'%010d'),'.glob.nc']);
    dicsurfdiag=rdmnc([modeldir,'/dic_surfDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
    dictave=rdmnc([modeldir,'/dic_tave.',num2str(start_iter,'%010d'),'.glob.nc']);
    
    ntsteps=size(tave.Ttave,4);
    
    % Prefer output from the diagnostics package over the timeave package
    if isfield(ocediag,'THETA')
        theta=ocediag.THETA; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        theta=tave.Ttave;
    end
    
    if isfield(ocediag,'SALT')
        salt=ocediag.SALT; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        salt=tave.Stave;
    end
    
    if isfield(dicflx,'TRAC01')
        dic=dicflx.TRAC01; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        dic=ptracers.dic;
    end
    
    if isfield(alkflx,'TRAC02')
        alk=alkflx.TRAC02; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        alk=ptracers.alk;
    end
    
    if isfield(po4flx,'TRAC03')
        po4=po4flx.TRAC03; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        po4=ptracers.po4;
    end
    
    if isfield(dopflx,'TRAC04')
        dop=dopflx.TRAC04; %.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        dop=ptracers.dop;
    end
    
    if isfield(ocediag,'UVEL')
        uvel=ocediag.UVEL.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        uvel=tave.uVeltave.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    end
    
    if isfield(ocediag,'VVEL')
        vvel=ocediag.VVEL.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        vvel=tave.vVeltave.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    end
    
    if isfield(ocediag,'WVEL')
        wvel=ocediag.WVEL.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    else
        wvel=tave.wVeltave.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    end
    
    if isfield(ocediag,'TOTTTEND')
        gtheta=ocediag.TOTTTEND./86400; % K/s not K/day
    else
        % Calculate tendency by the time differences
        gtheta(:,:,:,1:ntsteps-1)=diff(theta,1,4)/((60*60*24*360)/ntsteps);
        gtheta(:,:,:,ntsteps)=(theta(:,:,:,1)-theta(:,:,:,end))/((60*60*24*360)/ntsteps);
    end
    
    if isfield(ocediag,'TOTSTEND')
        gsalt =ocediag.TOTSTEND./86400; % psu/s not psu/day
    else
        % Calculate tendency by the time differences
        gsalt(:,:,:,1:ntsteps-1)=diff(salt,1,4)/((60*60*24*360)/ntsteps);
        gsalt(:,:,:,ntsteps)=(salt(:,:,:,1)-salt(:,:,:,end))/((60*60*24*360)/ntsteps);
    end
    
% For the passive tracers, physical tendency needs to be added to 
%   sources/sinks from gchem_forcing_sep  
    if isfield(dicflx,'gTr01') && isfield(dicflx,'DICGDIC')
        gdic=dicflx.gTr01+dicflx.DICGDIC; % Sum of PTRACER and DIC package diagnostics
    else
        % Calculate tendency by the time differences
        gdic(:,:,:,1:ntsteps-1)=diff(dic,1,4)/((60*60*24*360)/ntsteps);
        gdic(:,:,:,ntsteps)=(dic(:,:,:,1)-dic(:,:,:,end))/((60*60*24*360)/ntsteps);
    end
    
    if isfield(alkflx,'gTr02') && isfield(alkflx,'DICGALK')
        galk=alkflx.gTr02+alkflx.DICGALK; % Sum of PTRACER and DIC package diagnostics
    else
        % Calculate tendency by the time differences
        galk(:,:,:,1:ntsteps-1)=diff(alk,1,4)/((60*60*24*360)/ntsteps); 
        galk(:,:,:,ntsteps)=(alk(:,:,:,1)-alk(:,:,:,end))/((60*60*24*360)/ntsteps);
    end
    
    if isfield(po4flx,'gTr03') && isfield(po4flx,'DICGPO4')
        gpo4=po4flx.gTr03+po4flx.DICGPO4; % Sum of PTRACER and DIC package diagnostics
    else
        % Calculate tendency by the time differences
        gpo4(:,:,:,1:ntsteps-1)=diff(po4,1,4)/((60*60*24*360)/ntsteps); 
        gpo4(:,:,:,ntsteps)=(po4(:,:,:,1)-po4(:,:,:,end))/((60*60*24*360)/ntsteps);
    end

    % Load atmospheric CO2 data
    if isfield(dicsurfdiag,'DICATCO2') 
        atm_pco2=dicsurfdiag.DICATCO2.*1e6; % uatm
    elseif isfield(dicsurfdiag,'DICATCAR')
        atm_pco2=(dicsurfdiag.DICATCAR./1.77e20).*1e6; %uatm
    else
        atm_co2 = load([modeldir,'/dic_atmos.',num2str(start_iter,'%010d'),'.txt']);
        atm_co2(:,1)=(atm_co2(:,1)*43200)/(60*60*24*360); % timesteps
        atm_pco2=repmat(permute(atm_co2(2:end,3),[3,2,1]),[grid.nx,grid.ny,1]).*1e6; % uatm
    end
    
    gpco2(:,:,1:ntsteps-1)=diff(atm_pco2,1,3)/((60*60*24*360)/ntsteps); 
    gpco2(:,:,ntsteps)=(atm_pco2(:,:,1)-atm_pco2(:,:,end))/((60*60*24*360)/ntsteps);
    
%% Do Mixed layer depth calculations    
    if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
        mxldepth=nanmean(surfdiag.MXLDEPTH,4);
    else
        mxldepth=surfdiag.MXLDEPTH;
    end
    
    mld_array.maxmld=nanmax(squeeze(mxldepth).*repmat(grid.hfacc(:,:,1),[1,1,size(mxldepth,4)]),[],3);
    mld_array.minmld=nanmin(squeeze(mxldepth).*repmat(grid.hfacc(:,:,1),[1,1,size(mxldepth,4)]),[],3);
    mld_array.meanmld=nanmean(squeeze(mxldepth).*repmat(grid.hfacc(:,:,1),[1,1,size(mxldepth,4)]),3);
    mld_array.mld=squeeze(mxldepth).*repmat(grid.hfacc(:,:,1),[1,1,size(mxldepth,4)]);
    
    mld_array.kmaxmld=get_kave(cumsum(grid.dz),mld_array.maxmld);
    mld_array.kminmld=get_kave(cumsum(grid.dz),mld_array.minmld);
    mld_array.kmeanmld=get_kave(cumsum(grid.dz),mld_array.meanmld);
    mld_array.kmld=get_kave(cumsum(grid.dz),mld_array.mld);
    
    % set number of layers to integrate over
    if strcmp(kave_method,'arb_depth')
        %    (1=25m, 2=85m, 3=170m, 4=290m, 5=455m )        Coarse res model
        %    (1=10m, 2=20m, 3 = 35m, 4=55m, 5=75m 6=100m)   ECCO 1 degree model
        kave=get_kave(cumsum(grid.dz),kave);
    elseif strcmp(kave_method,'maxmld')
        %kave=get_kave(grid.zc,mld_array.maxmld);
        kave=mld_array.kmaxmld;
    elseif strcmp(kave_method,'minmld')
        %kave=get_kave(grid.zc,mld_array.minmld);
        kave=mld_array.kminmld;
    elseif strcmp(kave_method,'meanmld')
        %kave=get_kave(grid.zc,mld_array.meanmld);
        kave=mld_array.kmeanmld;
    elseif strcmp(kave_method,'mld')
        %kave=get_kave(grid.zc,mld_array.mld);
        kave=mld_array.kmld;
        disp('For seasonally varying MLD, need to include entrainment term, not implimented')
    else
        kave=ones(grid.nx,grid.ny);
    end
    
    if size(kave,3)~=size(theta,4)
        kave=repmat(kave,[1,1,size(theta,4)]);
    end

    % annual average if necessary
    if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
        % Needs a time dimension
        ntsteps=2;
        theta=repmat(nanmean(theta,4),[1,1,1,2]);
        salt=repmat(nanmean(salt,4),[1,1,1,2]);
        
        dic=repmat(nanmean(dic,4),[1,1,1,2]);
        alk=repmat(nanmean(alk,4),[1,1,1,2]);
        po4=repmat(nanmean(po4,4),[1,1,1,2]);
        dop=repmat(nanmean(dop,4),[1,1,1,2]);
        
        uvel=repmat(nanmean(uvel,4),[1,1,1,2]);
        vvel=repmat(nanmean(vvel,4),[1,1,1,2]);
        wvel=repmat(nanmean(wvel,4),[1,1,1,2]);
        
        gdic=repmat(nanmean(gdic,4),[1,1,1,2]);
        galk=repmat(nanmean(galk,4),[1,1,1,2]);
        gpo4=repmat(nanmean(gpo4,4),[1,1,1,2]);
    end
%%    
    if strcmp(obsid(1:min(13,end)),'mitgcm_fluxes')
        % Load advective and diffusive fluxes produced by the model.
        utflux=flxdiag.ADVx_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vtflux=flxdiag.ADVy_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wtflux=flxdiag.ADVr_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        utdiff=flxdiag.DFxE_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vtdiff=flxdiag.DFyE_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wtediff=flxdiag.DFrE_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wtidiff=flxdiag.DFrI_TH.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        
        usflux=flxdiag.ADVx_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vsflux=flxdiag.ADVy_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wsflux=flxdiag.ADVr_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        usdiff=flxdiag.DFxE_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vsdiff=flxdiag.DFyE_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wsediff=flxdiag.DFrE_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wsidiff=flxdiag.DFrI_SLT.*repmat(grid.hfacc,[1,1,1,ntsteps]);
                
        ucflux=dicflx.ADVxTr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vcflux=dicflx.ADVyTr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wcflux=dicflx.ADVrTr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        ucdiff=dicflx.DFxETr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vcdiff=dicflx.DFyETr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wcediff=dicflx.DFrETr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wcidiff=dicflx.DFrITr01.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        
        uaflux=alkflx.ADVxTr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vaflux=alkflx.ADVyTr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        waflux=alkflx.ADVrTr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        uadiff=alkflx.DFxETr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vadiff=alkflx.DFyETr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        waediff=alkflx.DFrETr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        waidiff=alkflx.DFrITr02.*repmat(grid.hfacc,[1,1,1,ntsteps]);

        upflux=po4flx.ADVxTr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vpflux=po4flx.ADVyTr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wpflux=po4flx.ADVrTr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        updiff=po4flx.DFxETr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        vpdiff=po4flx.DFyETr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wpediff=po4flx.DFrETr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        wpidiff=po4flx.DFrITr03.*repmat(grid.hfacc,[1,1,1,ntsteps]);
    elseif nc_isvar([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc'],'GM_U_RES') && ...
            nc_isvar([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc'],'GM_V_RES') && ...
            nc_isvar([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc'],'GM_W_RES')
        
        if ~isvar('timestep');
            timestep=mit_getparm([modeldir,'/data'],'deltaTClock');
            if isempty(timestep);
                timestep=43200;
            end
        end
        
        if ~isvar('advscheme');
            advscheme=mit_getparm([modeldir,'/data'],'tempAdvScheme');
            if isempty(advscheme);
                advscheme=2;
            end
        end
        
        if ~isvar('ptradvscheme')
            ptradvscheme=mit_getparm([modeldir,'/data.ptracers'],'PTRACERS_advScheme');
            if isempty(ptradvscheme)
                ptradvscheme=77;
            end
        end
        
        % annual average if necessary
        if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
            gmdiag=rdmnc([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
            ures=repmat(nanmean(gmdiag.GM_U_RES,4).*grid.hfacc,[1,1,1,2]);
            vres=repmat(nanmean(gmdiag.GM_V_RES,4).*grid.hfacc,[1,1,1,2]);
            wres=repmat(nanmean(gmdiag.GM_W_RES,4).*grid.hfacc,[1,1,1,2]);
        else
            gmdiag=rdmnc([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
            ures=gmdiag.GM_U_RES.*repmat(grid.hfacc,[1,1,1,ntsteps]);
            vres=gmdiag.GM_V_RES.*repmat(grid.hfacc,[1,1,1,ntsteps]);
            wres=gmdiag.GM_W_RES.*repmat(grid.hfacc,[1,1,1,ntsteps]);
        end
        
        %[uflux,vflux,wflux]=advection_flux(advscheme,timestep,dx,dy,dz,rax,ray,rac,...
        %                                  mask,uvel,vvel,wvel,scalar);
        [utflux,vtflux,wtflux]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(ures,'==',NaN,0),change(vres,'==',NaN,0),change(wres,'==',NaN,0),change(theta,'==',NaN,0));
        [usflux,vsflux,wsflux]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(ures,'==',NaN,0),change(vres,'==',NaN,0),change(wres,'==',NaN,0),change(salt,'==',NaN,0));
        [ucflux,vcflux,wcflux]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(ures,'==',NaN,0),change(vres,'==',NaN,0),change(wres,'==',NaN,0),change(dic,'==',NaN,0));
        [uaflux,vaflux,waflux]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(ures,'==',NaN,0),change(vres,'==',NaN,0),change(wres,'==',NaN,0),change(alk,'==',NaN,0));
        [upflux,vpflux,wpflux]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(ures,'==',NaN,0),change(vres,'==',NaN,0),change(wres,'==',NaN,0),change(po4,'==',NaN,0));
    else
        if ~isvar('timestep');
            timestep=mit_getparm([modeldir,'/data'],'deltaTClock');
            if isempty(timestep);
                timestep=43200;
            end
        end
        
        if ~isvar('advscheme');
            advscheme=mit_getparm([modeldir,'/data'],'tempAdvScheme');
            if isempty(advscheme);
                advscheme=2;
            end
        end
        
        if ~isvar('ptradvscheme')
            ptradvscheme=mit_getparm([modeldir,'/data.ptracers'],'PTRACERS_advScheme');
            if isempty(ptradvscheme)
                ptradvscheme=77;
            end
        end

        %[uflux,vflux,wflux]=advection_flux(advscheme,timestep,dx,dy,dz,rax,ray,rac,...
        %                                  mask,uvel,vvel,wvel,scalar);
        [utflux_eulmean,vtflux_eulmean,wtflux_eulmean]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(uvel,'==',NaN,0),change(vvel,'==',NaN,0),change(wvel,'==',NaN,0),change(theta,'==',NaN,0));
        [usflux_eulmean,vsflux_eulmean,wsflux_eulmean]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(uvel,'==',NaN,0),change(vvel,'==',NaN,0),change(wvel,'==',NaN,0),change(salt,'==',NaN,0));
        [ucflux_eulmean,vcflux_eulmean,wcflux_eulmean]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(uvel,'==',NaN,0),change(vvel,'==',NaN,0),change(wvel,'==',NaN,0),change(dic,'==',NaN,0));
        [uaflux_eulmean,vaflux_eulmean,waflux_eulmean]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(uvel,'==',NaN,0),change(vvel,'==',NaN,0),change(wvel,'==',NaN,0),change(alk,'==',NaN,0));
        [upflux_eulmean,vpflux_eulmean,wpflux_eulmean]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
            grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
            change(uvel,'==',NaN,0),change(vvel,'==',NaN,0),change(wvel,'==',NaN,0),change(po4,'==',NaN,0));
        
        % Whether or not using residual transport, need to load GM for diffusive fluxes
        gmdiag=rdmnc([modeldir,'/gmDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
        
        if strcmp(eddy_method,'residual')
            tave.ugmtrans=nan(grid.nx,grid.ny,grid.nz,ntsteps);
            tave.vgmtrans=nan(grid.nx,grid.ny,grid.nz,ntsteps);
            tave.wgmtrans=nan(grid.nx,grid.ny,grid.nz,ntsteps);
            
           if isfield(gmdiag,'GM_U_EDD') && isfield(gmdiag,'GM_V_EDD') && isfield(gmdiag,'GM_W_EDD')
               % annual average if necessary
               if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
                   uvel_eddy=repmat(nanmean(gmdiag.GM_U_EDD,4),[1,1,1,2]);
                   vvel_eddy=repmat(nanmean(gmdiag.GM_V_EDD,4),[1,1,1,2]);
                   wvel_eddy=repmat(nanmean(gmdiag.GM_W_EDD,4),[1,1,1,2]);
               else
                   uvel_eddy=gmdiag.GM_U_EDD;
                   vvel_eddy=gmdiag.GM_V_EDD;
                   wvel_eddy=gmdiag.GM_W_EDD;
               end
           else
               
               if strcmpi(mit_getparm([modeldir,'/data.gmredi'],'GM_AdvForm'),'true') % Using advective form of GM
                   % annual average if necessary
                   if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
                       gm_psix=repmat(nanmean(gmdiag.GM_PsiX,4),[1,1,1,2]);
                       gm_psiy=repmat(nanmean(gmdiag.GM_PsiY,4),[1,1,1,2]);
                   else
                       gm_psix=gmdiag.GM_PsiX;
                       gm_psiy=gmdiag.GM_PsiY;
                   end
                   
                   [uvel_eddy,vvel_eddy,wvel_eddy]=mit_gmvel(gm_psix,gm_psiy,...
                       grid.nx,grid.ny,grid.nz,ntsteps,grid.dxg,grid.dyg,grid.dz,grid.cmask,grid.umask,grid.vmask,grid.rac);
               else
                   % annual average if necessary
                   if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
                       gm_kwx=repmat(nanmean(gmdiag.GM_Kwx,4),[1,1,1,2]);
                       gm_kwy=repmat(nanmean(gmdiag.GM_Kwy,4),[1,1,1,2]);
                   else
                       gm_kwx=gmdiag.GM_Kwx;
                       gm_kwy=gmdiag.GM_Kwy;
                   end
                   [uvel_eddy,vvel_eddy,wvel_eddy]=mit_gmvel(gm_kwx,gm_kwy,...
                       grid.nx,grid.ny,grid.nz,ntsteps,grid.dxg,grid.dyg,grid.dz,grid.cmask,grid.umask,grid.vmask,grid.rac);
               end
           end
            
            [utflux_eddies,vtflux_eddies,wtflux_eddies]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
                grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
                change(uvel_eddy,'==',NaN,0),change(vvel_eddy,'==',NaN,0),change(wvel_eddy,'==',NaN,0),change(theta,'==',NaN,0));
            [usflux_eddies,vsflux_eddies,wsflux_eddies]=advection_flux(advscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
                grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
                change(uvel_eddy,'==',NaN,0),change(vvel_eddy,'==',NaN,0),change(wvel_eddy,'==',NaN,0),change(salt,'==',NaN,0));
            [ucflux_eddies,vcflux_eddies,wcflux_eddies]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
                grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
                change(uvel_eddy,'==',NaN,0),change(vvel_eddy,'==',NaN,0),change(wvel_eddy,'==',NaN,0),change(dic,'==',NaN,0));
            [uaflux_eddies,vaflux_eddies,waflux_eddies]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
                grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
                change(uvel_eddy,'==',NaN,0),change(vvel_eddy,'==',NaN,0),change(wvel_eddy,'==',NaN,0),change(alk,'==',NaN,0));
            [upflux_eddies,vpflux_eddies,wpflux_eddies]=advection_flux(ptradvscheme,timestep,grid.dxg,grid.dyg,grid.dz,...
                grid.rax,grid.ray,grid.rac,grid.umask,grid.vmask,grid.cmask,...
                change(uvel_eddy,'==',NaN,0),change(vvel_eddy,'==',NaN,0),change(wvel_eddy,'==',NaN,0),change(po4,'==',NaN,0));

            % Combine Eulerian-mean and eddy contributions
            utflux=utflux_eulmean+utflux_eddies; vtflux=vtflux_eulmean+vtflux_eddies; wtflux=wtflux_eulmean+wtflux_eddies;
            usflux=usflux_eulmean+usflux_eddies; vsflux=vsflux_eulmean+vsflux_eddies; wsflux=wsflux_eulmean+wsflux_eddies;
            ucflux=ucflux_eulmean+ucflux_eddies; vcflux=vcflux_eulmean+vcflux_eddies; wcflux=wcflux_eulmean+wcflux_eddies;
            upflux=upflux_eulmean+upflux_eddies; vpflux=vpflux_eulmean+vpflux_eddies; wpflux=wpflux_eulmean+wpflux_eddies;
            uaflux=uaflux_eulmean+uaflux_eddies; vaflux=vaflux_eulmean+vaflux_eddies; waflux=waflux_eulmean+waflux_eddies;
        else
            utflux=utflux_eulmean; vtflux=vtflux_eulmean; wtflux=wtflux_eulmean;
            usflux=usflux_eulmean; vsflux=vsflux_eulmean; wsflux=wsflux_eulmean;
            ucflux=ucflux_eulmean; vcflux=vcflux_eulmean; wcflux=wcflux_eulmean;
            upflux=upflux_eulmean; vpflux=vpflux_eulmean; wpflux=wpflux_eulmean;
            uaflux=uaflux_eulmean; vaflux=vaflux_eulmean; waflux=waflux_eulmean;                    
        end
        
        % Pre-allocate the "extra diagonals" with default values
        gm_k=mit_getparm([modeldir,'/data.gmredi'],'GM_isopycK');
        
        if isempty(gm_k);
            gm_k=mit_getparm([modeldir,'/data.gmredi'],'GM_background_K');
        end
        
        gm_kux=ones(size(theta)).*gm_k; % non-unity diagonal
        gm_kuz=zeros(size(theta));      % extra diagonal
        gm_kvy=ones(size(theta)).*gm_k; % non-unity diagonal
        gm_kvz=zeros(size(theta));      % extra diagonal
        
        if isvar('gmdiag.GM_Kvy'); gm_kvy=gmdiag.GM_Kvy.*repmat(grid.vmask,[1,1,1,size(gmdiag.GM_Kvy,4)]); end
        if isvar('gmdiag.GM_Kvz'); gm_kvz=gmdiag.GM_Kvz.*repmat(grid.vmask,[1,1,1,size(gmdiag.GM_Kvz,4)]); end
        if isvar('gmdiag.GM_Kux'); gm_kux=gmdiag.GM_Kux.*repmat(grid.umask,[1,1,1,size(gmdiag.GM_Kux,4)]); end
        if isvar('gmdiag.GM_Kuz'); gm_kuz=gmdiag.GM_Kuz.*repmat(grid.umask,[1,1,1,size(gmdiag.GM_Kuz,4)]); end
        gm_kwx=gmdiag.GM_Kwx.*repmat(grid.cmask,[1,1,1,size(gmdiag.GM_Kwx,4)]);
        gm_kwy=gmdiag.GM_Kwy.*repmat(grid.cmask,[1,1,1,size(gmdiag.GM_Kwy,4)]);
        gm_kwz=gmdiag.GM_Kwz.*repmat(grid.cmask,[1,1,1,size(gmdiag.GM_Kwz,4)]); % For diffusion calculations
        
        if mit_getparm([modeldir,'/data.pkg'],'useKPP')
            disp('loading KPP diffusion values')
            kppdiag=rdmnc([modeldir,'/kppDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
            diffkz=(kppdiag.KPPdiffS)+repmat(mit_getparm([modeldir,'/data'],'diffKzS'),size(kppdiag.KPPdiffS))+gm_kwz;
        elseif mit_getparm([modeldir,'/data'],'ivdc_kappa')>0
            convadj=ocediag.CONVADJ.*repmat(grid.hfacc,[1,1,1,size(ocediag.CONVADJ,4)]);
            diffkz=(convadj.*mit_getparm([modeldir,'/data'],'ivdc_kappa'))+repmat(mit_getparm([modeldir,'/data'],'diffKzS'),size(convadj))+gm_kwz;
        else
            diffkz=repmat(mit_getparm([modeldir,'/data'],'diffKzS'),size(tave.Ttave))+gm_kwz;
        end
        
        if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
            gm_kux=repmat(nanmean(gm_kux,4),[1,1,1,2]);
            gm_kuz=repmat(nanmean(gm_kuz,4),[1,1,1,2]);
            gm_kvy=repmat(nanmean(gm_kvy,4),[1,1,1,2]);
            gm_kvz=repmat(nanmean(gm_kvz,4),[1,1,1,2]);
            gm_kwx=repmat(nanmean(gm_kwx,4),[1,1,1,2]);
            gm_kwy=repmat(nanmean(gm_kwy,4),[1,1,1,2]);
            gm_kwz=repmat(nanmean(gm_kwz,4),[1,1,1,2]);
            diffkz=repmat(nanmean(diffkz,4),[1,1,1,2]);
        end
    end % End of online/offline fluxes distinction
    
%% Surface Fluxes of Heat, Salt, CO2, virtual fluxes of DIC and ALK
    heatflux=squeeze(surfdiag.TFLUX).*repmat(grid.hfacc(:,:,1),[1,1,size(surfdiag.TFLUX,4)]); % W m-2
    theta_surf_flux=heatflux./(rhoconst*cp); % K m s-1
    
    salt_surf_flux=squeeze(surfdiag.SFLUX).*repmat(grid.hfacc(:,:,1),[1,1,size(surfdiag.SFLUX,4)])./rhoconst; % g m-2 s-1 * m3 kg-1 => Psu m s-1
   
    dic_tot_flux=squeeze(dicsurfdiag.DICTFLX).*repmat(grid.hfacc(:,:,1),[1,1,size(dicsurfdiag.DICTFLX,4)]).*grid.dz(1); % mol m-2 s-1
    dic_co2_flux=squeeze(dicsurfdiag.DICCFLX).*repmat(grid.hfacc(:,:,1),[1,1,size(dicsurfdiag.DICCFLX,4)]); % mol m-2 s-1
    
    if isfield(dicsurfdiag,'DICVCFLX')
        dic_vrt_flux=squeeze(dicsurfdiag.DICVCFLX).*repmat(grid.hfacc(:,:,1),[1,1,size(dicsurfdiag.DICVCFLX,4)]); % mol m-2 s-1
    else
        area=repmat(grid.rac,[1,1,12]);
        surf_salt = squeeze(salt(:,:,1,:)); smean=nansum(surf_salt(:).*area(:))./nansum(area(:));
        surf_dic  = squeeze(dic(:,:,1,:));  cmean=nansum(surf_dic (:).*area(:))./nansum(area(:));
        dic_vrt_flux=salt_surf_flux.*(cmean/smean);
    end
    
    if isfield(dicsurfdiag,'DICAFLX')
        alk_vrt_flux=squeeze(dicsurfdiag.DICAFLX).*repmat(grid.hfacc(:,:,1),[1,1,size(dicsurfdiag.DICAFLX,4)]).*grid.dz(1); % mol m-2 s-1
    else
        area=repmat(grid.rac,[1,1,12]);
        surf_salt = squeeze(salt(:,:,1,:)); smean=nansum(surf_salt(:).*area(:))./nansum(area(:));
        surf_alk  = squeeze(alk(:,:,1,:));  amean=nansum(surf_alk (:).*area(:))./nansum(area(:));
        alk_vrt_flux=salt_surf_flux.*(amean/smean);
    end    

    % Load KPP non-local fluxes and add to the advection terms
    if strcmpi(mit_getparm([modeldir,'/data.pkg'],'useKPP'),'true') && ...
            nc_isvar([modeldir,'/kppDiag.',num2str(start_iter,'%010d'),'.glob.nc'],'KPPg_TH');
        disp('Loading KPP Non-local vertical transports')
        kppdiag=rdmnc([modeldir,'/kppDiag.',num2str(start_iter,'%010d'),'.glob.nc']);
        
       if strcmp(obsid(1:min(13,end)),'mitgcm_fluxes')
            wtflux=wtflux+kppdiag.KPPg_TH.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPg_TH,4)]);
            wsflux=wsflux+kppdiag.KPPg_SLT.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPg_SLT,4)]);
%             wcflux=wcflux+kppdiag.KPPgTr01.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr01,4)]);
%             waflux=waflux+kppdiag.KPPgTr02.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr02,4)]);
%             wpflux=wpflux+kppdiag.KPPgTr03.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr03,4)]);
        else
            % Calculate offline the non-local fluxes: -rA.*(KPPdiffKzT/S.*KPPghat).*surf_flux
            kppg_th =-repmat(grid.rac,[1,1,grid.nz,size(theta_surf_flux,3)]).*kppdiag.KPPghatK.*repmat(permute(theta_surf_flux,[1,2,4,3]),[1,1,grid.nz,1]);
            kppg_slt=-repmat(grid.rac,[1,1,grid.nz,size(salt_surf_flux ,3)]).*kppdiag.KPPghatK.*repmat(permute(salt_surf_flux ,[1,2,4,3]),[1,1,grid.nz,1]);
            
            if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
                wtflux=wtflux+repmat(nanmean(kppg_th .*repmat(grid.cmask,[1,1,1,size(kppg_th ,4)]),4),[1,1,1,2]);
                wsflux=wsflux+repmat(nanmean(kppg_slt.*repmat(grid.cmask,[1,1,1,size(kppg_slt,4)]),4),[1,1,1,2]);
                %        wcflux=wcflux+nanmean(kppdiag.KPPgTr01.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr01,4)]),4);
                %        waflux=waflux+nanmean(kppdiag.KPPgTr02.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr02,4)]),4);
                %        wpflux=wpflux+nanmean(kppdiag.KPPgTr03.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr03,4)]),4);
            else
                wtflux=wtflux+kppg_th .*repmat(grid.cmask,[1,1,1,size(kppg_th ,4)]);
                wsflux=wsflux+kppg_slt.*repmat(grid.cmask,[1,1,1,size(kppg_slt,4)]);
                %        wcflux=wcflux+kppdiag.KPPgTr01.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr01,4)]);
                %        waflux=waflux+kppdiag.KPPgTr02.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr02,4)]);
                %        wpflux=wpflux+kppdiag.KPPgTr03.*repmat(grid.cmask,[1,1,1,size(kppdiag.KPPgTr03,4)]);
            end
        end
    else
        disp('No non-local transports from KPP')
    end

%% Biological production terms
    dic_bio_ave=dictave.dic_BIO_ave.*repmat(grid.hfacc,[1,1,1,size(dictave.dic_BIO_ave,4)]);
    dic_pflux_ave=dictave.dic_pflux_ave.*repmat(grid.hfacc,[1,1,1,size(dictave.dic_pflux_ave,4)]);
    dic_car_ave=dictave.dic_CAR_ave.*repmat(grid.hfacc,[1,1,1,size(dictave.dic_CAR_ave,4)]);
    
    if strcmp(obsid(1:min(13,end)),'mitgcm_annual')
        theta_surf_flux=repmat(nanmean(theta_surf_flux,3),[1,1,2]);
        salt_surf_flux =repmat(nanmean(salt_surf_flux,3),[1,1,2]);
        dic_tot_flux   =repmat(nanmean(dic_tot_flux,3),[1,1,2]);
        dic_co2_flux   =repmat(nanmean(dic_co2_flux,3),[1,1,2]);
        dic_vrt_flux   =repmat(nanmean(dic_vrt_flux,3),[1,1,2]);
        alk_vrt_flux   =repmat(nanmean(alk_vrt_flux,3),[1,1,2]);
        
        dic_bio_ave    =repmat(nanmean(dic_bio_ave,4),[1,1,1,2]);
        dic_pflux_ave  =repmat(nanmean(dic_pflux_ave,4),[1,1,1,2]);
        dic_car_ave    =repmat(nanmean(dic_car_ave,4),[1,1,1,2]);
        
        gtheta         =repmat(nanmean(gtheta,4),[1,1,1,2]);
        gsalt          =repmat(nanmean(gsalt,4),[1,1,1,2]);
    end

% calculate mixed layer integrated tendencies
[dthetadt]=vertically_integrate(grid.dz,kave,grid.cmask,gtheta);
[dsaltdt ]=vertically_integrate(grid.dz,kave,grid.cmask,gsalt );
[ddicdt  ]=vertically_integrate(grid.dz,kave,grid.cmask,gdic  );
[dalkdt  ]=vertically_integrate(grid.dz,kave,grid.cmask,galk  );
[dpo4dt  ]=vertically_integrate(grid.dz,kave,grid.cmask,gpo4  );
   
[po4_bio_flux]=vertically_integrate(grid.dz,kave,grid.cmask,-dic_bio_ave+dic_pflux_ave+max(KDOPRemin.*dop,0));
[dic_carb_flux]=vertically_integrate(grid.dz,kave,grid.cmask,dic_car_ave);

% biological activity
dic_bio_flux=Rcp.*po4_bio_flux; % mol m-2 s-1
% carbonate
alk_carb_flux=2.*dic_carb_flux; % mol m-2 s-1
% biological activity
alk_bio_flux=-Rnp.*po4_bio_flux; % mol m-2 s-1

%clear dictave ptracers
%% Calculate component transport, divergence, losses and gains
% Find advective fluxes using v.c (e.g. from data) or import
% the correct model fields.

% Calculate divergence of advective fluxes (calculated offline or online) and integrate vertically by kave levels
[theta_adv_horz,theta_adv_vert]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
    change(utflux,'==',NaN,0),change(vtflux,'==',NaN,0),change(wtflux,'==',NaN,0));
[salt_adv_horz ,salt_adv_vert ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
    change(usflux,'==',NaN,0),change(vsflux,'==',NaN,0),change(wsflux,'==',NaN,0));
[dic_adv_horz  ,dic_adv_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
    change(ucflux,'==',NaN,0),change(vcflux,'==',NaN,0),change(wcflux,'==',NaN,0));
[alk_adv_horz  ,alk_adv_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
    change(uaflux,'==',NaN,0),change(vaflux,'==',NaN,0),change(waflux,'==',NaN,0));
[po4_adv_horz  ,po4_adv_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
    change(upflux,'==',NaN,0),change(vpflux,'==',NaN,0),change(wpflux,'==',NaN,0));

% Calculate diffusive fluxes and then their divergence
switch lower(vertmix)
    case('')
        % Calculate divergence of diffusive fluxes from the model
        [theta_diff_horz,theta_diff_vert]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(utdiff,'==',NaN,0),change(vtdiff,'==',NaN,0),change(wtediff+wtidiff,'==',NaN,0));
        [salt_diff_horz,salt_diff_vert ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(usdiff,'==',NaN,0),change(vsdiff,'==',NaN,0),change(wsediff+wsidiff,'==',NaN,0));
        [dic_diff_horz,dic_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(ucdiff,'==',NaN,0),change(vcdiff,'==',NaN,0),change(wcediff+wcidiff,'==',NaN,0));
        [alk_diff_horz,alk_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(uadiff,'==',NaN,0),change(vadiff,'==',NaN,0),change(waediff+waidiff,'==',NaN,0));
        [po4_diff_horz,po4_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(updiff,'==',NaN,0),change(vpdiff,'==',NaN,0),change(wpediff+wpidiff,'==',NaN,0));
    case('diffusion')
        % Get the explicit vertical diffusion
        [utdiff,vtdiff,wtediff]=diffusion_flux(grid.dxc,grid.dyc,grid.dzc,...
            grid.umask,grid.vmask,grid.cmask,grid.rax,grid.ray,grid.rac,...
            gm_kux,gm_kuz,gm_kvy,gm_kvz,gm_kwx,gm_kwy,change(theta,'==',NaN,0));
        [usdiff,vsdiff,wsediff]=diffusion_flux(grid.dxc,grid.dyc,grid.dzc,...
            grid.umask,grid.vmask,grid.cmask,grid.rax,grid.ray,grid.rac,...
            gm_kux,gm_kuz,gm_kvy,gm_kvz,gm_kwx,gm_kwy,change(salt,'==',NaN,0));
        [ucdiff,vcdiff,wcediff]=diffusion_flux(grid.dxc,grid.dyc,grid.dzc,...
            grid.umask,grid.vmask,grid.cmask,grid.rax,grid.ray,grid.rac,...
            gm_kux,gm_kuz,gm_kvy,gm_kvz,gm_kwx,gm_kwy,change(dic,'==',NaN,0));
        [uadiff,vadiff,waediff]=diffusion_flux(grid.dxc,grid.dyc,grid.dzc,...
            grid.umask,grid.vmask,grid.cmask,grid.rax,grid.ray,grid.rac,...
            gm_kux,gm_kuz,gm_kvy,gm_kvz,gm_kwx,gm_kwy,change(alk,'==',NaN,0));
        [updiff,vpdiff,wpediff]=diffusion_flux(grid.dxc,grid.dyc,grid.dzc,...
            grid.umask,grid.vmask,grid.cmask,grid.rax,grid.ray,grid.rac,...
            gm_kux,gm_kuz,gm_kvy,gm_kvz,gm_kwx,gm_kwy,change(po4,'==',NaN,0));
        
        % Calculate implicit vertical mixing using  MEX
        % This code iterates a few times (that's fine according to Large et al
        % 1994) and actually outputs the fluxes at each iteration, so you can check
        % to see if when they converge
        timestep=43200;
        iteration=12; % iteration of the implicit matrix inversion to return
        wtidiff=implicit_diffusion_flux(timestep,iteration,grid.rac,grid.dz,grid.dzc,grid.cmask,...
            diffkz,theta);
        %        change(diffkz,'==',NaN,0),change(theta,'==',NaN,0));
        wsidiff=implicit_diffusion_flux(timestep,iteration,grid.rac,grid.dz,grid.dzc,grid.cmask,...
            diffkz,salt);
        %        change(diffkz,'==',NaN,0),change(salt,'==',NaN,0));
        wcidiff=implicit_diffusion_flux(timestep,iteration,grid.rac,grid.dz,grid.dzc,grid.cmask,...
            diffkz,dic);
        %        change(diffkz,'==',NaN,0),change(dic,'==',NaN,0));
        waidiff=implicit_diffusion_flux(timestep,iteration,grid.rac,grid.dz,grid.dzc,grid.cmask,...
            diffkz,alk);
        %        change(diffkz,'==',NaN,0),change(alk,'==',NaN,0));
        wpidiff=implicit_diffusion_flux(timestep,iteration,grid.rac,grid.dz,grid.dzc,grid.cmask,...
            diffkz,po4);
        %        change(diffkz,'==',NaN,0),change(po4,'==',NaN,0));
        
        % Calculate divergence of diffusive fluxes
        [theta_diff_horz,theta_diff_vert]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(utdiff,'==',NaN,0),change(vtdiff,'==',NaN,0),change(wtediff+wtidiff,'==',NaN,0));
        [salt_diff_horz,salt_diff_vert ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(usdiff,'==',NaN,0),change(vsdiff,'==',NaN,0),change(wsediff+wsidiff,'==',NaN,0));
        [dic_diff_horz,dic_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(ucdiff,'==',NaN,0),change(vcdiff,'==',NaN,0),change(wcediff+wcidiff,'==',NaN,0));
        [alk_diff_horz,alk_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(uadiff,'==',NaN,0),change(vadiff,'==',NaN,0),change(waediff+waidiff,'==',NaN,0));
        [po4_diff_horz,po4_diff_vert  ]=calc_divergence(grid.dz,grid.volc,kave,grid.cmask,...
            change(updiff,'==',NaN,0),change(vpdiff,'==',NaN,0),change(wpediff+wpidiff,'==',NaN,0)); %#ok<*ASGLU>
    otherwise
        theta_diff_vert=theta_adv_horz.*0;
        salt_diff_vert=theta_adv_horz.*0;
        dic_diff_vert=theta_adv_horz.*0;
        alk_diff_vert=theta_adv_horz.*0;
        po4_diff_vert=theta_adv_horz.*0;
        
        theta_diff_horz=theta_diff_vert.*0;
        salt_diff_horz=theta_diff_vert.*0;
        dic_diff_horz=theta_diff_vert.*0;
        alk_diff_horz=theta_diff_vert.*0;
        po4_diff_horz=theta_diff_vert.*0;
end

%% Plot budget terms and combine to produce CO2 flux decomposition
%
% Csat is a linear function of T so dCsat ~ dDIC/dTheta*dTheta
% dDIC/dT = -8.976 x 10-6  mol kg-1 K-1 from Mick, Similar to Goodwin and Lenton (2007) who use 0.01 mol/m3/K
%ddicdt = -8.5069e-6 ; % Calculated from co2sys routine @ S=34.828, Alk=2304.4, pco2=280
%dDICdS = -5.7435e-6 ; % Calculated from co2sys routine @ T=18.792, Alk=2304.4, pco2=280
%dDICdALK = 0.80526 ;  % Calculated from co2sys routine @ T=18.792, S=34.828, pco2=280
% Slightly different depending on atm CO2...

% Convert to umol/kg for calculations (which produce mol kg-1 per x) then convert back to mol/m3
[ddicdth,ddicds,ddicdalk,ddicdpco2]=calc_ddicdvar(grid,theta(:,:,1,:),salt(:,:,1,:),alk(:,:,1,:).*molperm3_2_umolperkg,po4(:,:,1,:).*molperm3_2_umolperkg,atm_pco2);
ddicdth=ddicdth.*rhoconst;
ddicds=ddicds.*rhoconst;
%ddicdalk=1 ;ddicdalk; % Just to see what happens really.
ddicdpco2(:,2)=ddicdpco2(:,2).*rhoconst;

% I have outsourced this to keep consistent
co2_flux_steadystate

%% Save output to a mat file
fieldnms = who('alk','alk_*','dic','dic_*','po4','po4_*','salt','salt_*','theta','theta_*',... % Temperature budget components
               '*dt','flux_*','biop','adv*','diff*','ugrad*','cres*','csat*','cflux','atm_pco2',... % reconstruction components
               'KDOPRemin','Rcaco3','Rcn','Rcp','Rno','Rnp','Rpo','Rsip','cp','grid','kave','mld_array','ntsteps',...                              
               'rhoconst','molperm3_2_umolperkg','perl_2_perkg','perm3_2_perl','seconds_per_year','umolperkgmpers_2_molperm2peryr',...
               'eddy_method','ent_method','kave_method','modeldir','obsid','vertmix');
                         
for ii=1:length(fieldnms)
    tmps.(fieldnms{ii})=eval([fieldnms{ii}]);
end

eval([obsid,'=tmps;']); clear tmps
fname=[obsid,'_diags.mat'];
disp(['saving ',obsid,' to ',fname])
save(fname,obsid)

