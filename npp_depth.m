function float_npp_depth(dpath, fpath, movtspath, npath, float_ids)
%%-------------------------------------------------------------------------
%
% Filename:				get_npp_from_argo.m
%
% Author:				Jacqueline S. Long
% Created:				10-Dec-2020
%
% Changed on:               06-23-2022  by N.V. Buzby   
% SHORT CHANGE DESCRIPTION: changed PAR data to be from MODIS data (rather
%                           than from OSU site), matched to float
%
%--------------------------------------------------------------------------
%
% Description:
%		This program reads in bio-optical float data and calculates net
%		primary production through depth using the Carbon-based Production
%		Model (CbPM, Behrenfeld et al., XXXX, Westberry 2008).
%
% INPUT VARIABLES:
%   dpath: depth averaged data path
%   fpath: float data path - original float data 
%   movtspath: Path to satellite data extracted at float location
%   npath: Path to save NPP float data structure
%
% Input parameters from float files:
% All parameters below should be 1xN, where N is the full profile length
%		- bbp			Backscatter [m^-1]
%		- Chl			Chlorophyll [mg m^-3]
%       - mld			Depth of the Mixed Layer [m]
%       - PAR			An estimate of photosynthetically available radiation at the surface
%       - posvec		[Latitude(ºN) Longitude(ºW)];
%       - dt			Datetime in Matlab SDN format
%       - press			Pressure or depth [dbar, m], from zmax to surface
%
% OPTIONAL INPUTS:
%   float_ids:  array of floats of interest, if not included will compute
%               for all existing downloaded float data
% OUTPUT:
%   will save calculated NPP as a structure, labeld with float ID
%   1-m binned data from 1 to 200 meters
%		- npp
%		- mu
%		- Chl:C
%		- irr
%
% Functions called:
% AustinPetzold_CBpM2
% get_doy
% floatCbPM_JSLV3_6
%
%% --------------------------------------------------------------------------
% Define input chl and bbp fields
bbp_fn1 = 'bbp700';
bbp_fn2 = bbp_fn1;
chl_fn = 'chla';

% Load data
%are float ids specified?
if nargin > 4
else %if not create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
        float_ids{i} = od_files_char(1:7);
    end
end

%% --------------------------    Program start   --------------------------
%                 No modifications should be made below this line
%--------------------------------------------------------------------------

zvec = [1:1:200];
null = 1;

%% Loop through floats
for i = 1:length(float_ids)
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = char(float_ids(i));
    else
        floatID = num2str(float_ids(i));
    end  
    load([fpath floatID '.mat'],'f')
   
    % Load Data
    load([dpath floatID '_depavg_data.mat']) %Load MLD data, mldmean struct
    load([fpath floatID '.mat']) %Load full profile data, f struct
    
    % Get original fieldnames to use later in interpolation loop
    fn = fieldnames(f);
    [prof pidx] = unique(f.profile);
    % Call 'f' something else so that we can save as 'f'
    fl = f; clear f
    % create variable to list profiles for which can't calc NPP
    bad_profs = 0;
    bad_var = {NaN,NaN};
    %% Loop through Profiles
    for p = 1:length(prof)
        % Get index of current profile
        current_idx = find(fl.profile == prof(p));
        
        if isempty(find(~isnan(fl.press(current_idx)), 1))
            warning('%s has empty pressure values for profile %s',...
                floatID, num2str(prof(p)))
        else
            %% PAR / IRR
            load([movtspath floatID '.mat']) %using float-matched satellite PAR values
            clear irr; irr = uf.par(p);
            clear uf
            
            %% Position Vector
            [~, moidx, ~, doy] = get_doy(odmean.date(p));
            clear posvec; posvec = [odmean.lat(p), odmean.lon(p)];

            %% MLD
            mld = mldmean.mld(p);

            %% OD
            od = odmean.od(p);

            %% Chl
            % Extrapolate the shallowest value to get a surface estimate
            clear CHLsurf_idx; CHLsurf_idx = find(~isnan(fl.(chl_fn)(current_idx)),1);
            if isempty(CHLsurf_idx) %then the whole profile is NaN
                Chl = NaN;
            else
                % Add it to chl vector
                clear tmp; tmp = [fl.(chl_fn)(current_idx(CHLsurf_idx));fl.(chl_fn)(current_idx)];
                nandex = ~isnan(fl.press(current_idx));
                tmpP = [fl.press(current_idx(nandex)); 1];
                % Sample points must be unique and ascending...
                [tmpP,uidx] = unique(tmpP);
                clear tmpidx; tmpidx = find(~isnan(tmp(uidx)));
                if length(tmpP(tmpidx)) < 2 || length(tmp(uidx(tmpidx))) < 2
                    Chl = NaN;
                else
                    Chl = interp1(tmpP(tmpidx), tmp(uidx(tmpidx)),zvec);
                    Chl(Chl < 0) = 0;
                end
            end
            %% Bbp
            % Extrapolate the shallowest value to get a surface estimate
            clear BBPsurf_idx; BBPsurf_idx = find(~isnan(fl.(bbp_fn1)(current_idx)),1);
            if isempty(BBPsurf_idx) %then the whole profile is NaN
                bbp = NaN;
            elseif isnan(nansum(fl.(bbp_fn2)(current_idx)))
                bbp = NaN;
            else
                clear BBPsurf_val; BBPsurf_val = fl.(bbp_fn1)(current_idx(BBPsurf_idx));
                % Add it to bbp vector
                clear tmp; tmp = [fl.(bbp_fn2)(current_idx);BBPsurf_val];
                if ~exist('uidx','var') %if empty Chl profile, need to establish unique/ascending order
                    nandex = ~isnan(fl.press(current_idx));
                    tmpP = [fl.press(current_idx(nandex)); 1];
                    [tmpP,uidx] = unique(tmpP);
                end
                clear tmpidx; tmpidx = find(~isnan(tmp(uidx))); %uidx refers to pressure values from Chl calcs
                if length(tmpidx) < 2
                    bbp = NaN;
                else
                    bbp700 = interp1(tmpP(tmpidx), tmp(uidx(tmpidx)),zvec); %tmpP also refers to previous press values
                    %Calculate bbp 470 from 700 > if using Graff relationship >
                    %this needs to also be modified in the calculationg of NPP
                    %function below
                    %             bbp = bbp700.*(470/700).^-1;  %Assume lambda^-1 dependence of scattering
                    %Calculate bbp 443 from 700 > if using CbPM (Behrenfeld) relationship
                    bbp = bbp700.*(443/700).^-1;  %Assume lambda^-1 dependence of scattering
                end
            end
            %% Kd
            % Using OD-averaged Chl
            if isnan(Chl)
                k490=NaN;
            else
                k490 = 0.0166 + 0.0773.*(nanmean(Chl(find(zvec < od)))).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
                k490(k490<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold, L&W Table 3.16
            end

            %% CALCULATE NPP

            % To save as original float 'f' structure size
            pp = p;
            % Add surface value to pressure
            new_press = [fl.press(current_idx)];
            % Need to add a data point for the index
            new_idx = [current_idx(1)+(p - 1):1:current_idx(end)+p];

            if (mld > 0) && nansum(Chl)>0 && ~isnan(irr) && ~isnan(k490)...
                    && nansum(bbp)>0 && ~isnan(sum(posvec))
                
                % Estimate NPP using CbPM code:
                clear out; [out] = floatCbPM_JSLV3_6(mld,Chl,bbp,k490,irr,posvec,doy,zvec,od);
                % Log parameters to save to float structure
                % Make changes here if you want to change how the output is saved (e.g., if
                % you want it to be save like the float data is
                f.npp_vec(:,pp) = interp1(zvec,out.ppz,new_press);
                f.mu(:,pp) = interp1(zvec,out.mu,new_press);
                f.chlz(:,pp) = interp1(zvec,Chl,new_press);
                f.PhytC(:,pp) = interp1(zvec,out.carbon,new_press);
                f.ChlC(:,pp) = interp1(zvec,out.chl_C,new_press);
                f.Igrid(:,pp) = interp1(zvec,out.parz,new_press);
                f.prcnt(:,pp) = interp1(zvec,out.prcnt,new_press);
                f.inpp.intNPP(p,1) = out.inpp;
                
                %Sensitivity Analysis
                perturb = {'Chl','bbp','label';1.5,1,'A';2,1,'B';4,1,'C';1,1.5,'A';1,2,'B';1,4,'C'};
                for b=2:length(perturb)
                        clear out; [out] = floatCbPM_JSLV3_6(mld,Chl*perturb{b,1},bbp*perturb{b,2},k490,irr,posvec,doy,zvec,od);
                        if b < 5
                            f.inpp.perturb.chl.(perturb{b,3})(p,1)= out.inpp;
                        else
                            f.inpp.perturb.bbp.(perturb{b,3})(p,1)=out.inpp;
                    end
                end

            else %MLD is 0 so not calculating NPP
                %organize into array to assign bad variable to profile
                if ~nansum(Chl)>0
                    bad_var(p,:) = {prof(p),'Chl'};
                elseif isnan(irr)
                    bad_var(p,:) = {prof(p),'irr'};
                elseif isnan(k490)
                    bad_var(p,:) = {prof(p),'kd490'};
                elseif ~nansum(bbp)>0
                    bad_var(p,:) = {prof(p),'bbp'};
                elseif isnan(sum(posvec))
                    bad_var(p,:) = {prof(p),'posvec'};
                elseif isnan(mld)
                    if ~nansum(fl.sal(:,p))
                        bad_var(p,:) = {prof(p),'sal'};
                    elseif ~nansum(fl.temp(:,p))
                        bad_var(p,:) = {prof(p),'temp'};
                    end
                else
                    bad_var(p,:) = {prof(p),'something else'};
                end   
                disp(['one of the input vars is NaN - ',char(bad_var(p,2))]);
                null = 0;
            end
            if null == 0 %then set to NaN for this profile
                rep_len = length(fl.date);
                f.prof_check(:,pp) = NaN*ones(rep_len,1);
                f.npp_vec(:,pp) = NaN*ones(rep_len,1);
                f.mu(:,pp) = NaN*ones(rep_len,1);
                f.chlz(:,pp) = NaN*ones(rep_len,1);
                f.PhytC(:,pp) = NaN*ones(rep_len,1);
                f.ChlC(:,pp) = NaN*ones(rep_len,1);
                f.Igrid(:,pp) = NaN*ones(rep_len,1);
                f.prcnt(:,pp) = NaN*ones(rep_len,1);
                f.inpp.intNPP(p,1) = NaN;
                
                %add NaN into perturbation vectors as well
                f.inpp.perturb.chl.A(p,1) = NaN;
                f.inpp.perturb.chl.B(p,1) = NaN;
                f.inpp.perturb.chl.C(p,1) = NaN;
                f.inpp.perturb.bbp.A(p,1) = NaN;
                f.inpp.perturb.bbp.B(p,1) = NaN;
                f.inpp.perturb.bbp.C(p,1) = NaN;
            else
            end
            f.date_vec(:,pp) = odmean.date(p);
            f.inpp.date(p,1) = odmean.date(p);
            f.inpp.mld(p,1) = mld;
            f.inpp.od(p,1) = odmean.od(p);
            f.inpp.zeu(p,1) = zeumean.zeu(p);
            % Loop through the original float fields to make the profiles line
            % up - this will just add a NaN value for alot these variables
            % at the surface but could change to "surface most"
            for k = 1:length(fn)
                current_var = char(fn(k));
                if isnan(nansum(fl.(current_var)(current_idx)))
                    f.(current_var)(:,pp) = NaN;
                else
                    if strcmp(current_var,'profile') || strcmp(current_var,'date') || strcmp(current_var,'lon') || strcmp(current_var,'lat')
                        f.(current_var)(:,pp) = fl.(current_var)(current_idx(1));
                    elseif strcmp(current_var,'press')
                        f.(current_var)(:,pp) = new_press;
                    else
                         % Sample points must be unique and ascending...
                        nandex = ~isnan(fl.press(current_idx));
                        clear tmpP uidx; [tmpP,uidx] = unique(fl.press(current_idx(nandex)));
                        clear tmp; tmp = fl.(current_var)(current_idx);
                        clear tmpidx; tmpidx = find(~isnan(tmp(uidx)));
                        if length(tmpidx) < 2 % If this profile of data doesn't have enough ~isnan values
                            f.(current_var)(:,pp) = nan(length(new_press),1);
                        else
                            f.(current_var)(:,pp) = interp1(tmpP(tmpidx),tmp(uidx(tmpidx)),new_press);
                        end
                    end
                end
            end

            %check_struct(f, {})
            
            clearvars -except pp floatID reg f fn fl bad_var mldmean zeumean odmean float_ids...
                uf prof bad_profs float_files od_files_char col bpath fpath movtspath zvec ...
                spath ppath bbp_fn1 bbp_fn2 chl_fn plist bad_var npath dpath
            null = 1; %reset null
        end
    end
    % Find latest float data to log
    f.info = {['Float ID: ' floatID]; ...
        ['Data processed on ' datestr(now)]; ...
        ['Float data time range = ' datestr(min(odmean.date)) ' to ' datestr(max(odmean.date))]};

    % Save npp data
    save([npath floatID '_CbPM.mat'],'f');
end

