function npp_surface(pp, dpath, fpath, movtspath, npath, platform, float_ids)
%Program to calculate NPP using the satellite algorithm (CbPM)
%with input data:
% Chl and bbp (from satellite obs - 8 days averaged)
% Nitrate (for zno3) from WOA2018 data - 1deg climatology bins extracted at
% sat lat/lon
% MLD from HYCOM data (monthly climatology data)
%
% Written July 24, 2020 Jacki Long (MBARI)
%
% Needs OD-averaged float data from: "B_extract_data_at_depth_horizons.m"
%
% INPUT VARIABLES:
%   pp: building off path
%   dpath: base path - depth averaged data 
%   fpath: float data path - original float data 
%   movtspath: Path to satellite data extracted at float location
%   npath: Path to save NPP float data structure
%   platform: 'float' or 'sat', based on what platform want to use for npp
%   calculation
%
% OPTIONAL INPUTS:
%   float_ids:  array of floats of interest, if not included will compute
%               for all existing downloaded float data
% Output:
%   will save calculated NPP as a structure 'snpp', labeld with float ID
%   1-m binned data from 1 to 200 meters
%		- inpp / npp_vec
%		- mu
%		- Chl:C
%		- irr
%       - PhytC
% Functions called:
% AustinPetzold_CBpM2
% get_doy
% floatCbPM_JSLV3_2_movTSdata

% Need to define input chl and bbp fields
bbp_fn1 = 'bp700';
bbp_fn2 = bbp_fn1;
chl_fn = 'chla';

%% Load data
% Load data
%are float ids specified?
if nargin > 6
else %if not create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
        float_ids{i} = od_files_char(1:7);
    end
end

%Load the MLD data from HYCOM (3D struct "M" on monthly data)
load([pp '/input/input_for_sat_comparison/HYCOM_MLD_071422.mat'])

%Load the zno3 data from WOA analysis (3D struct "Z" on monthly clims)
load([pp '/input/input_for_sat_comparison/zno3_monthly_clim.mat'])
Z = zno3; clear zno3

%%  Loop through all floats and each profile
null = 1;
for i = 1:length(float_ids)
    % What float are we on?
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = float_ids{i};
    else
        floatID = num2str(float_ids(i));
    end  
    f = load([fpath floatID '.mat']); %Load full profile data ("f" struct)
    load([dpath floatID '_depavg_data.mat']); %Load MLD data ("mldmean" and "mldstd" structs)
    uf = load([movtspath floatID '.mat']); %Load movTS satellite data ("uf" struct)
    
    [prof profidx] = unique(f.profile);
    %% Loop through Profiles
    for p = 1:length(uf.profile)
        % Current idx of full profile data, if needed:
        current_idx = find(f.profile == prof(p));
        % Get PAR data (same for both float & satellite estimate)
        clear irr; irr = uf.par(p);
        
        % Get general info of this profile from the float obs data
        plat = f.lat(p);
        plon = f.lon(p);
        pdate = f.date(p);
        
        % Find daylength at this latitude and DOY
        if strcmp(platform, 'sat') % For satellite data:
            [~, moidx, ~, doy] = get_doy(pdate);
            daylength = day_length(doy,plat);

        elseif strcmp(platform, 'float') % For float data:
            [~, moidx, ~, doy] = get_doy(odmean.date(p));
            plat = odmean.lat(p); %serving posvec purpose
            daylength = day_length(doy,plat);
        end
        %% MLD
        % 7/25/22 changed to using float MLD b/c HYCOM dataset ends in 2020
        % MLD found using Holte & Talley method
        mld = mldmean.mld(p);

%         if strcmp(platform, 'float') %For float data:
%         % This is Holte & Talley method (also westberry option?)
%             mld = mldmean.mld(p);
%         
%         elseif strcmp(platform, 'sat') % For satellite data:
%         % Get MLD estimate from HYCOM data ('M') & Interpolate to lat lon and date of this profile
%           mld = interp3(M.lon, M.lat, M.date, M.mld, plon, plat, pdate); %DOESN'T WORK WHEN FLOAT DATE > HYCOM (see HYCOM_MLD_dataupdate.mat)
%         end
        %% OD
        od = odmean.od(p);
        zvec = [1:1:round(od)]';
        %% Zno3
        if strcmp(platform, 'float') %For float data:
            % Zno3 depth is defined as the depth where nitrate + nitrite exceeded 0.5 uM
            % 0.5 umol/L * rho => umol/kg
            % gsw_rho returns kg/m => /1000 => L
            clear frho; frho = sw_dens(f.sal(current_idx),f.temp(current_idx),f.press(current_idx))./1000;
            clear tmpno3; tmpno3 = (f.no3(current_idx) ./ frho); %Now in uM
            % In this region, I think with this definition of the
            % nitracline, we are always below it (e.g., surface values are
            % always >0.5uM
            clear tmpidx; tmpidx = find(tmpno3 >= 0.5,1,'last');
            if ~isempty(tmpidx)
                clear zno3; zno3 = find(f.press(current_idx(tmpidx)));
            else %use previous zno3 profile value as estimate
                disp(['No Zno3 found at profile ' num2str(p) ' float ' floatID])
            end
        
        elseif strcmp(platform, 'sat') %For satellite data:
            % Using WOA monthly clim data ('Z') Extract zno3 
            % (calculated from the depth resolved grid data in a previous step calc_zno3_WOA.m)
             zno3 = interp2(Z.lon,Z.lat,Z.data(:,:,moidx), plon, plat);
        end
        %% Chl
        if strcmp(platform, 'float') %For float data:
            %Chl = odmean.est_chla(p); 
            
            %Changing to extrapolating Chl to surface, interpolating to 1-m bins and averaging...
            %Surface-most not NaN value,
            clear CHLsurf_idx; CHLsurf_idx = find(~isnan(f.(chl_fn)(current_idx)),1,'last');
            if isempty(CHLsurf_idx) %then the whole profile is NaN
                Chl = NaN;
            else
                % Add it to chl vector
                clear tmp; tmp = [f.(chl_fn)(current_idx);f.(chl_fn)(current_idx(CHLsurf_idx))];
                tmpP = [f.press(current_idx); 1];
                % Interpolate to 1-m bins
                Chl_OD = interp1(tmpP, tmp,zvec);
                Chl_OD(Chl_OD < 0) = 0;
                Chl = nanmean(Chl_OD);
            end
        
        elseif strcmp(platform, 'sat') %For satellite data:
            Chl = uf.chlor_a(p);
        end
        %% Bbp
        if strcmp(platform, 'float') % For float data
            %bbp700 = odmean.bbp_tot(p);
        
            % Changing to extrapolating Bbp to surface, interpolating to 1-m bins and averaging...
            clear BBPsurf_idx; BBPsurf_idx = find(~isnan(f.(bbp_fn1)(current_idx)),1,'last');
            if isempty(BBPsurf_idx) %then the whole profile is NaN
                bbp = NaN;
            elseif isnan(nansum(f.(bbp_fn2)(current_idx)))
                bbp = NaN;
            else
                clear BBPsurf_val; BBPsurf_val = f.(bbp_fn1)(current_idx(BBPsurf_idx));
                % Add it to chl vector
                clear tmp; tmp = [f.(bbp_fn2)(current_idx);BBPsurf_val];
                tmpidx = find(~isnan(tmp));
                bbp700_OD = interp1(tmpP(tmpidx), tmp(tmpidx),zvec);
                %Calculate bbp 443 from 700 > if using CbPM relationship (also Graff option)
                bbp_OD = bbp700_OD.*(443/700).^-1;  %Assume lambda^-1 dependence of scattering
                bbp = nanmean(bbp_OD);
            end
        
        elseif strcmp(platform, 'sat') %For satellite data:
            bbp = uf.bbp_443_gsm(p);
        end
        %% Kd
        if strcmp(platform, 'float') %For float data:
        %compute Kd490, KdPAR, and propagate satellite PAR through the surface
        %Assuming this gives kd490 in units 1/m, for CbPM code requirement
            k490 = 0.0166 + 0.0773.*(Chl).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
            k490(k490<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold, L&W Table 3.16

        elseif strcmp(platform, 'sat') %For satellite data:
         k490 = uf.Kd_490(p);
        end
        %% Calculate NPP
        if (mld > 0) && ~isnan(Chl) && ~isnan(irr)  && ~isnan(k490) && ~isnan(bbp)
            
            clear out; [out] = floatCbPM_JSLV3_2_movTSdata(mld,Chl,bbp,k490,zno3,irr,daylength);

            % Log parameters to save to float structure
            snpp.depth(:,p) = out.z;
            snpp.inpp.intNPP(p) = out.inpp;
            snpp.npp_vec(:,p) = out.ppz;
            snpp.mu(:,p) = out.mu;
            snpp.chlz(:,p) = out.chlz;
            snpp.PhytC(:,p) = out.carbon;
            snpp.ChlC(:,p) = out.chl_C;
            snpp.Igrid(:,p) = out.parz;
            snpp.prcnt(:,p) = out.prcnt;
            snpp.mld([1:200],p) = mld;
            snpp.od([1:200],p) = odmean.od(p);
            snpp.zeu([1:200],p) = zeumean.zeu(p);
            snpp.lat([1:200],p) = mldmean.lat(p);
            snpp.lon([1:200],p) = mldmean.lon(p);
            
            %Sensitivity Analysis
%             perturb = {'Chl','bbp','label';1.5,1,'A';2,1,'B';4,1,'C';1,1.5,'A';1,2,'B';1,4,'C'};
%             for b=2:length(perturb)
%                 clear out; [out] = floatCbPM_JSLV3_2_movTSdata(mld,Chl*perturb{b,1},bbp*perturb{b,2},k490,zno3,irr,daylength);
%                 if b < 5
%                 	snpp.inpp.perturb.chl.(perturb{b,3})(p)= out.inpp;
%                 else
%                 	snpp.inpp.perturb.bbp.(perturb{b,3})(p)=out.inpp;
%                 end
%             end
            
        else %MLD is 0 so not calculating NPP
            disp([floatID ' profile ' num2str(p) ', MLD = 0 or ODchl = NaN']);
            null = 0;
        end
        
        if null == 0 %then set to NaN for this profile
            snpp.depth(:,p) = NaN.*ones(200,1);
            snpp.inpp.intNPP(p) = NaN;
            snpp.npp_vec(:,p) = NaN.*ones(200,1);
            snpp.mu(:,p) = NaN.*ones(200,1);
            snpp.chlz(:,p) = NaN.*ones(200,1);
            snpp.PhytC(:,p) = NaN;%.*ones(200,1);
            snpp.ChlC(:,p) = NaN.*ones(200,1);
            snpp.Igrid(:,p) = NaN.*ones(200,1);
            snpp.prcnt(:,p) = NaN.*ones(200,1);
            snpp.mld([1:200],p) = NaN;
            snpp.od([1:200],p) = NaN;
            snpp.zeu([1:200],p) = NaN;
            snpp.lat([1:200],p) = NaN;
            snpp.lon([1:200],p) = NaN;
        else
        end
        snpp.date([1:200],p) = odmean.date(p);
        null = 1; %reset null
    end
    %% Save npp data
    isfolder(npath);
    save([npath floatID '_' platform '_surface_CbPM.mat'],'-struct','snpp');
    
    clear snpp mldmean mldstd floatID od_files_char odmean odstd profile_num
end