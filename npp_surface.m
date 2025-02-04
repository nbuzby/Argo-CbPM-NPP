function npp_surface(pp, dpath, fpath, movtspath, npath, platform, float_ids, matchup_type, chl_correction)
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
%   float_ids:      array of floats of interest, if not included will compute
%                   for all existing downloaded float data
%   matchup_types:  specifies matchup type for satellite vaules (closest
%                   point to float profile, or mean. default is mean.
%   chl_correction: determine whether to use original float chlorophyll
%                   values or satellite-corrected values (Y/N). if not
%                   included will use corrected values
%
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

%set matchup type if not specified, default is closest
if nargin < 8
    matchup_type = 'med';
    chl_correction = 'Y';
end

% Need to define input chl fields
if strcmp(chl_correction, 'Y')
    chl_fn = 'chl_corr';
else
    chl_fn = 'chla';
end
null = 1;


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
for i = 1:length(float_ids)
    % What float are we on?
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = float_ids{i};
    else
        floatID = num2str(float_ids(i));
    end  
    f = load([fpath floatID '.mat']); %Load full profile data ("f" struct)
    load([dpath floatID '_depavg_data.mat'],'mldmean','odmean','zeumean'); %Load depth-avged data 
    if strcmp(platform, 'sat')
        uf = load([movtspath floatID '.mat']); %Load movTS satellite data ("uf" struct)
    end
    [prof, ~] = unique(f.profile);
    %% Loop through Profiles
    for p = 1:width(f.profile)
        if strcmp(platform, 'sat') % For satellite data:
            %Determine closest point for satellite estimates
            clear closest_idx;
            [~,closest_idx] = min(distance(uf.lat(p),uf.lon(p),uf.match.lat(:,p),...
            uf.match.lon(:,p)));
            pdate = uf.match.date(p);
             % Lat/Lon pair & Find daylength at this latitude and DOY 
            if strcmp(matchup_type,'closest')
                plat = uf.match.lat(closest_idx,p);
                plon =  uf.match.lon(closest_idx,p); 
            elseif strcmp(matchup_type,'mean')
                plat = uf.lat(p);
                plon = uf.lon(p);
            end
            [~, moidx, ~, doy] = get_doy(pdate);
            daylength = day_length(doy,plat);
            

        elseif strcmp(platform, 'float') % For float data:
            plat = f.lat(p);
            plon = f.lon(p);
            pdate = f.date(p);
            [~, moidx, ~, doy] = get_doy(odmean.date(p));
            daylength = day_length(doy,plat);
        end

        %% PAR
        if strcmp(platform, 'float') %For float data:
            clear irr;
            % if isfield(f,'irr_int')
            %     if p < length(odmean.irr_int)
            %         irr = odmean.irr_int(p);
            %     else 
            %         irr = NaN;
            %     end
            % else
                if length(f.par)>1 && p<length(f.par)+1
                    irr = f.par(p);
                else
                    irr = NaN;
                end
            %end

        elseif strcmp(platform, 'sat') % For satellite data
            if strcmp(matchup_type,'closest')
                clear irr; irr= uf.par(closest_idx,p);
            elseif strcmp(matchup_type,'mean')
                clear irr; irr= mean(uf.par(:,p),'omitnan');
            end
        end
        %% MLD
        % 7/25/22 changed to using float MLD b/c HYCOM dataset ends in 2020
        % MLD found using de Boyer & Montegut method
        mld = f.mld(1,p);
   
%         strcmp(platform, 'sat') % For satellite data:
%         % Get MLD estimate from HYCOM data ('M') & Interpolate to lat lon and date of this profile
%           mld = interp3(M.lon, M.lat, M.date, M.mld, plon, plat, pdate); %DOESN'T WORK WHEN FLOAT DATE > HYCOM (see HYCOM_MLD_dataupdate.mat)
        %% OD
        od = odmean.od(p);
        zvec = (1:1:round(od))';
        %% Zno3
        if strcmp(platform, 'float') %For float data:
            if isfield(f,'no3') && ~isnan(sum(f.no3(:,p)))
                % Zno3 depth is defined as the depth where nitrate + nitrite exceeded 0.5 uM
                % 0.5 umol/L * rho => umol/kg
                % gsw_rho returns kg/m => /1000 => L
                clear frho; frho = sw_dens(f.sal(:,p),f.temp(:,p),f.press(:,p))./1000;
                clear tmpno3; tmpno3 = (f.no3(:,p) ./ frho); %Now in uM
                % In this region, I think with this definition of the
                % nitracline, we are always below it (e.g., surface values are
                % always >0.5uM
                clear tmpidx; tmpidx = find(tmpno3 >= 0.5,1,'last');
                if ~isempty(tmpidx)
                    clear zno3; zno3 = find(f.press(tmpidx,p));
                else %use previous zno3 profile value as estimate
                    disp(['No Zno3 found at profile ' num2str(p) ' float ' floatID])
                end
            else
                zno3 = interp2(Z.lon,Z.lat,Z.data(:,:,moidx), plon, plat);
            end
        
        elseif strcmp(platform, 'sat') %For satellite data:
            % Using WOA monthly clim data ('Z') Extract zno3 
            % (calculated from the depth resolved grid data in a previous step calc_zno3_WOA.m)
             zno3 = interp2(Z.lon,Z.lat,Z.data(:,:,moidx), plon, plat);
        end
        %% Chl
        if strcmp(platform, 'float') %For float data:
            if strcmp(chl_correction,'Y')
                Chl = odmean.(chl_fn).(matchup_type)(p); 
            else
                Chl = odmean.(chl_fn)(p);
            end
        elseif strcmp(platform, 'sat') %For satellite data:
            if strcmp(matchup_type,'closest')
                clear Chl; Chl= uf.chlor_a(closest_idx,p);
            elseif strcmp(matchup_type,'mean')
                clear Chl; Chl= mean(uf.chlor_a(:,p),'omitnan');
            end
        end
        %% Bbp
        if strcmp(platform, 'float') % For float data
            bbp = odmean.bbp700(p)*(443/700)^-1;
        
        elseif strcmp(platform, 'sat') %For satellite data:
            if strcmp(matchup_type,'closest')
                clear bbp; bbp= uf.bbp_443(closest_idx,p);
            elseif strcmp(matchup_type,'mean')
                clear bbp; bbp= mean(uf.bbp_443(:,p),'omitnan');
            end
        end
        %% Kd
        if strcmp(platform, 'float') %For float data:
        %compute Kd490, KdPAR, and propagate satellite PAR through the surface
        %Assuming this gives kd490 in units 1/m, for CbPM code requirement
            k490 = 0.0166 + 0.0773.*(Chl).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
            k490(k490<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold, L&W Table 3.16

        elseif strcmp(platform, 'sat') %For satellite data:
            k490 = uf.Kd_490(closest_idx,p);
            if strcmp(matchup_type,'closest')
                clear k490; k490= uf.Kd_490(closest_idx,p);
            elseif strcmp(matchup_type,'mean')
                clear k490; k490= mean(uf.Kd_490(:,p),'omitnan');
            end
        end
        %% Calculate NPP
        if (mld > 1) && ~isnan(Chl) && ~isnan(irr)  && ~isnan(k490) && ~isnan(bbp)
            
            clear out; [out] = floatCbPM_JSLV3_2_movTSdata(mld,Chl,bbp,k490,zno3,irr,daylength);

            % Log parameters to save to float structure
            snpp.depth(:,p) = out.z;
            snpp.inpp.intNPP(p) = out.inpp;
            snpp.npp_vec(:,p) = out.ppz;
            snpp.mu(:,p) = out.mu;
            snpp.chlz(:,p) = out.chlz;
            snpp.PhytC(:,p) = out.carbon;
            snpp.ChlC(:,p) = out.chl_C;
            snpp.parz(:,p) = out.parz;
            snpp.Igrid(:,p) = out.parz;
            snpp.prcnt(:,p) = out.prcnt;
            snpp.mld(1:200,p) = mld;
            snpp.od(1:200,p) = odmean.od(p);
            snpp.zeu(1:200,p) = zeumean.zeu(p);
            snpp.lat(1:200,p) = mldmean.lat(p);
            snpp.lon(1:200,p) = mldmean.lon(p);
            
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
            %organize into array to assign bad variable to profile
            if ~nansum(Chl)>0
                snpp.bad_var(p,:) = {prof(p),'Chl'};
            elseif isnan(irr)
                snpp.bad_var(p,:) = {prof(p),'irr'};
            elseif isnan(k490)
                snpp.bad_var(p,:) = {prof(p),'kd490'};
            elseif ~nansum(bbp)>0
                snpp.bad_var(p,:) = {prof(p),'bbp'};
            elseif isnan(mld)
                    snpp.bad_var(p,:) = {prof(p),'mld'};
            else
                snpp.bad_var(p,:) = {prof(p),'something else'};
            end   
            %disp(['For float ' floatID ' profile ' num2str(prof(p)) ' one of the input vars is NaN - ',char(snpp.(matchup_type).bad_var(p,2))]);
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
            snpp.parz(:,p) = NaN.*ones(200,1);
            snpp.Igrid(:,p) = NaN.*ones(200,1);
            snpp.prcnt(:,p) = NaN.*ones(200,1);
            snpp.mld(1:200,p) = NaN;
            snpp.od(1:200,p) = NaN;
            snpp.zeu(1:200,p) = NaN;
            snpp.lat(1:200,p) = NaN;
            snpp.lon(1:200,p) = NaN;
        else
        end
        snpp.date(1:200,p) = odmean.date(1,p);
        null = 1; %reset null
    end
        snpp.info = {['Float ID: ' floatID]; ...
        ['Data processed on ' datestr(now)]; ...
        ['Float data time range = ' datestr(min(min(odmean.date))) ' to ' datestr(max(max(odmean.date)))];...
        ['Matchup method= ' matchup_type]};

    %% Save npp data
    % isfolder(npath);
    % if isfile([npath floatID '_' platform '_surface_CbPM.mat'])
    %     save([npath floatID '_' platform '_surface_CbPM.mat'],'-struct','snpp','-append');
    % else
    if strcmp(chl_correction,'Y') 
        tmp = platform;
    else
        tmp = [platform '_uncorr'];
    end

    % if isfield(f,'irr_int')
    %     save([npath floatID '_' tmp '_irr_surface_CbPM.mat'],'-struct','snpp');
    % else
        save([npath floatID '_' tmp '_surface_CbPM.mat'],'-struct','snpp');
   % end
    %end
    
    clear snpp mldmean mldstd floatID od_files_char odmean odstd profile_num
end