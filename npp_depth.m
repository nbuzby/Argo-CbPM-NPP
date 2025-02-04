function npp_depth(dpath, fpath, npath, float_ids, chl_correction,par_data)
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
%   float_ids:      array of floats of interest, if not included will compute
%                   for all existing downloaded float data
%   chl_correction: determine whether to use original float chlorophyll
%                   values or satellite-corrected values (Y/N). if not
%                   included will use corrected values
%   par_data:       choose what PAR data source to use: satellite matchup
%                   (sat), depth resolved ML product from Herve's group
%                   (prod), or float radiometer data (float).
%                   if not included will use satellite matchup
%
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
% Load data
%are float ids specified? and should corrected chlorophyll be used?

if nargin < 6
    par_data = 'sat';
elseif nargin < 5
    chl_correction = 'Y';
    par_data = 'sat';
elseif nargin < 4
    %if not create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
        float_ids{i} = od_files_char(1:7);
    end
    chl_correction = 'Y';
    par_data = 'sat';
end

% Define input chl field
if chl_correction == 'Y'
    chl_fn = 'chla_corr';
else
    chl_fn = 'chla';
end

%% --------------------------    Program start   --------------------------
%                 No modifications should be made below this line
%--------------------------------------------------------------------------

zvec = 1:1:200;
null = 1;
perturb = {'Chl','bbp','par','label';1.5,1,1,1;2,1,1,2;4,1,1,3;1,1.5,1,1;1,2,1,2;1,4,1,3;1,1,1.5,1;1,1,2,2;1,1,4,3};

%% Loop through floats
for i = 1:length(float_ids)
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = char(float_ids(i));
    else
        floatID = num2str(float_ids(i));
    end  
   
    % Load Data
    load([dpath floatID '_depavg_data.mat'],'odmean','zeumean') %Load dep-avged data (mld, od, zeu structs)
    f = load([fpath floatID '.mat']); %Load full profile data, f struct
    
    prof = 1:width(f.profile);
    % Call 'f' something else so that we can save as 'f'
    fl = f; clear f

    %% Loop through Profiles
    for p = 1:length(prof)
        if isempty(find(~isnan(fl.press(:,p)), 1))
            warning('%s has empty pressure values for profile %s',...
                floatID, num2str(prof(p)))
            f.date_vec(:,p) = odmean.date(p);
            f.inpp.date(p,1) = odmean.date(p);
            f.inpp.mld(p,1) = mld;
            f.inpp.od(p,1) = odmean.od(p);
            f.inpp.zeu(p,1) = zeumean.zeu(p);
            f.bad_var(p,:) = {prof(p),'press'};
        else
            %% Define Input Variables for CbPM
            [~, ~, ~, doy] = get_doy(odmean.date(p));
            clear posvec; posvec = [odmean.lat(p), odmean.lon(p)]; %Position Vector
            
            if strcmp(par_data,'sat') && length(fl.par)>1
                clear irr; irr = fl.par(p); %using float-matched satellite PAR values
            elseif strcmp(par_data,'prod')
                clear irr; irr = fl.PAR_prod(1:200,p); %using ML PAR product from Renosh et al., 2023 (DOI: 10.3390/rs15245663)
                irr(irr<0)=NaN;
            elseif strcmp(par_data,'float')
                clear irr; irr = fl.irr_int(1:200,p);
                irr(irr<0)=NaN;
            else
                irr = NaN;
            end

            mld = fl.mld(1,p); 
            od = odmean.od(p); 
            Chl = fl.(chl_fn)(1:200,p); Chl(isnan(Chl)) = 0; %can't have NaNs at surface
            bbp = fl.bbp700(1:200,p).*(443/700).^-1; bbp(isnan(bbp)) = 0; %Calculate bbp 443 from 700 > if using CbPM (Behrenfeld) relationship (Assume lambda^-1 dependence of scattering)
        
            %Kd (Using OD-averaged Chl)
            if isnan(Chl)
                k490=NaN;
            else
                k490 = 0.0166 + 0.0773.*(nanmean(Chl(zvec < od))).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
                k490(k490<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold, L&W Table 3.16
            end

            %% CALCULATE NPP
            if (mld > 0) && nansum(Chl)>0 && nansum(irr)>0 && ~isnan(k490)...
                    && nansum(bbp)>0 && ~isnan(sum(posvec))
                
                % Estimate NPP using CbPM code:
                if strcmp(par_data,'sat')
                    clear out; [out] = floatCbPM_JSLV3_6(mld,Chl,bbp,k490,irr,posvec,doy,zvec,od);
                elseif strcmp(par_data,'prod') || strcmp(par_data,'float')
                    clear out; [out] = floatCbPM_JSLV3_6_PARZ(Chl,bbp,irr,posvec,doy,zvec,od);
                end
                % Log parameters to save to float structure
                f.npp_vec(:,p) = out.ppz;
                f.parz(:,p) = out.parz;
                f.mu(:,p) = out.mu;
                f.chlz(:,p) = Chl;
                f.PhytC(:,p) = out.carbon;
                f.ChlC(:,p) = out.chl_C;
                f.Igrid(:,p) = out.parz;
                f.prcnt(:,p) = out.prcnt;
                f.inpp.intNPP(p,1) = out.inpp;
                f.bad_var(p,:) = {prof(p),'None'};
                
                %Sensitivity Analysis
                for b=2:length(perturb)
                    if strcmp(par_data,'sat')
                        clear out; [out] = floatCbPM_JSLV3_6(mld,Chl*perturb{b,1},bbp*perturb{b,2},k490,irr*perturb{b,3},posvec,doy,zvec,od);
                    elseif strcmp(par_data,'prod')
                        clear out; [out] = floatCbPM_JSLV3_6_PARZ(Chl*perturb{b,1},bbp*perturb{b,2},irr*perturb{b,3},posvec,doy,zvec,od);
                    end
                    
                    out.inpp(out.inpp==0) = NaN;
                    if b < 5
                        f.inpp.perturb.chl(p,perturb{b,4})= out.inpp;
                    elseif b>4 && b<8
                        f.inpp.perturb.bbp(p,perturb{b,4})=out.inpp;
                    else
                        f.inpp.perturb.par(p,perturb{b,4})=out.inpp;
                    end
                end

            else %MLD is 0 so not calculating NPP
                %organize into array to assign bad variable to profile
                if ~nansum(Chl)>0
                    f.bad_var(p,:) = {prof(p),'Chl'};
                elseif isnan(irr)
                    f.bad_var(p,:) = {prof(p),['irr_' par_data]};
                elseif isnan(k490)
                    f.bad_var(p,:) = {prof(p),'kd490'};
                elseif ~nansum(bbp)>0
                    f.bad_var(p,:) = {prof(p),'bbp'};
                elseif isnan(sum(posvec))
                    f.bad_var(p,:) = {prof(p),'posvec'};
                elseif isnan(mld)
                    if ~nansum(fl.sal(:,p))
                        f.bad_var(p,:) = {prof(p),'sal'};
                    elseif ~nansum(fl.temp(:,p))
                        f.bad_var(p,:) = {prof(p),'temp'};
                    else
                        f.bad_var(p,:) = {prof(p),'MLD'};
                    end
                else
                    f.bad_var(p,:) = {prof(p),'something else'};
                end
                %disp(['For float ' floatID ' profile ' num2str(prof(p)) ' one of the input vars is NaN - ',char(f.bad_var(p,2))]);
                null = 0;
            end
            if null == 0 %then set to NaN for this profile
                rep_len = length(zvec);
                f.prof_check(:,p) = NaN*ones(rep_len,1);
                f.npp_vec(:,p) = NaN*ones(rep_len,1);
                f.parz(:,p) = NaN*ones(rep_len,1);
                f.mu(:,p) = NaN*ones(rep_len,1);
                f.chlz(:,p) = NaN*ones(rep_len,1);
                f.PhytC(:,p) = NaN*ones(rep_len,1);
                f.ChlC(:,p) = NaN*ones(rep_len,1);
                f.Igrid(:,p) = NaN*ones(rep_len,1);
                f.prcnt(:,p) = NaN*ones(rep_len,1);
                f.inpp.intNPP(p,1) = NaN;
                
                %add NaN into perturbation vectors as well
                f.inpp.perturb.chl(p,:) = NaN;
                f.inpp.perturb.bbp(p,:) = NaN;
                f.inpp.perturb.par(p,:) = NaN;
            end
        end
        f.date_vec(p) = odmean.date(1,p);
        f.inpp.date(p,1) = odmean.date(1,p);
        f.inpp.mld(p,1) = mld;
        f.inpp.od(p,1) = odmean.od(p);
        f.inpp.zeu(p,1) = zeumean.zeu(p);
        f.press(:,p) = zvec;
        f.chl(:,p) = Chl;
        f.bbp(:,p) = bbp;

        
        clearvars -except floatID f fn fl bad_var mldmean zeumean odmean float_ids...
            uf prof bad_profs float_files od_files_char col bpath fpath movtspath zvec ...
            spath ppath bbp_fn1 bbp_fn2 chl_fn bad_var npath dpath perturb chl_correction par_data
        null = 1; %reset null
    end
    % Find latest float data to log
    f.info = {['Float ID: ' floatID]; ...
        ['Data processed on ' datestr(now)]; ...
        ['Float data time range = ' datestr(min(min(odmean.date))) ' to ' datestr(max(max(odmean.date)))]};

    % Save npp data
    if chl_correction == 'Y'
        chl_corr = '';
    elseif chl_correction == 'N'
        chl_corr = '_uncorrected';
    end
    if strcmp(par_data,'prod')
        par_name = '_parproduct';
    elseif strcmp(par_data,'sat')
        par_name = '';
    elseif strcmp(par_data,'float')
        par_name = 'float_rad';
    end
    filename = ['_CbPM' chl_corr par_name];
    
    if isfile(filename)
        save([npath floatID filename],'-struct','f','-append'); 
    else 
        save([npath floatID filename],'-struct','f');
    end
end

