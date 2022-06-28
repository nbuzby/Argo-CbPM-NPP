function [mld, mldmean, mldstd,desc] = calc_mld_float(f, method, mean_calcs, ftype)
%
% [mld, mldmean, mldstd] = calc_mld_float(f, method, mean_calcs)
%
% Find MLD with option to average data within MLD
%
% INPUT VARIABLES:
%   f: float struture
%   method: 'BY' or 'HT' (see below)
%   mean_calcs: boolean where
%       1 = calculate MLD averages
%       0 = do not calculate MLD averages
%
% OUTPUT VARIABLES:
%   mld: mixed layer depth
%   mldmean: structure of MLD averaged variables
%   mldstd: structure of MLD averaged variable StDevs
%
% de Boyer Montégut et al. (2004)
%   'BY' method
%       0.03 kg/m3 dnesity increase relative to 10m
%       0.2 C temperature decrease relative to 10m
%
% Holte & Talley (2009) - https://doi.org/10.1175/2009JTECHO543.1
%   'HT' method
%       Cite as:
%           Holte, J., L. D. Talley, J. Gilson, and D. Roemmich (2017),
%           An Argo mixed layer climatology and database, Geophys. Res. Lett., 44,
%           5618?5626, doi:10.1002/2017GL073426."
%
% FUNCTIONS CALLED: 
%	HT method: findmld.m downloaded from http://mixedlayer.ucsd.edu/
%
% FUNCTION HISTORY:
%	Written by Jacki Long (MBARI) Oct 2019
% 	Updated by Andrea Fassbender April 2020
%   May 12 2020 (JL), removed error causing data overlap, also changed to search 
%   for entire character string to match variables using strmatch rather than contains- J. Long
%   Last Update:
%   July 24 2020 (JL), Added desc output for units. Now Zeu and OD averaged
%   data are also saved as structs all in same file under
%   "B_extract_data_at_depth_horizons.m"
%% Define Variables

% Predefine struct variables
mldmean = struct;
mldstd = struct;

fieldstr  = {'temp','sal','sigt','oxy','oxysat','no3','chla','est_chla','qc_chla','LPCTD_adj_chla','bp700',...
             'bbp_tot','poc','ph_insitu','ta_liar','dic_liar','co2_liar','cdom','sat_par','dens'};
fieldstr2 = {'date','profile','lon','lat'};
             
units = {'Deg C','PSU','kg m^-3','umol kg^-1','%','umol kg^-1','ug L^-1','ug L^-1','ug L^-1','ug L^-1','m^-1',...
         'm^-1','mmol m^-3','total scale','umol kg^-1','umol kg^-1','uatm','ppb','E m^-2 d^-1','rho'};
     
% Find unique profiles
[profile,~] = unique(f.profile);

% Pre-allocate 
mld = NaN(length(profile),1);

% Calculate density
if ~isfield(f,'dens')
    f.dens =  sw_dens(f.sal,f.temp,f.press); %kg/m3
else
end

%% Calculate MLD & MLD Averages if Specified

for i = 1:length(profile)
    
    % Isolate non NaN profile data
    clear p currentidx salNaN tempNaN presNaN noNaN; 
    p = find(f.profile==profile(i));
    salNaN = ~isnan(f.sal(p));
    tempNaN = ~isnan(f.temp(p));
    presNaN = ~isnan(f.press(p));
    currentidx = p(find(salNaN.*tempNaN.*presNaN >0));
    
    % Define input parameters for findmld.m
    if ftype == 'Argo' %data is already in (pres) ascending order
        sal = f.sal(currentidx);
        temp = f.temp(currentidx);
        pres = f.press(currentidx);
        dens = f.dens(currentidx);
    else %have to flip data
        sal = flipud(f.sal(currentidx));
        temp = flipud(f.temp(currentidx));
        pres = flipud(f.press(currentidx));
        dens = flipud(f.dens(currentidx));
    end
    
    % Find profiles that are all NaN values
    if sum([nansum(sal),nansum(temp),nansum(pres),nansum(dens)])==0
    else
        % Isolate non NaN profile data
        clear p currentidx salNaN tempNaN presNaN noNaN; 
        p = find(f.profile==profile(i));
        salNaN = ~isnan(f.sal(p));
        tempNaN = ~isnan(f.temp(p));
        presNaN = ~isnan(f.press(p));
        currentidx = p(find(salNaN.*tempNaN.*presNaN >0));
        % Calculate MLD
        if strcmp(method,'BY') == 1
%             % de Boyer Montégut et al. method:
%             
%             % Find 10m temperature and 0.2 temperature change
%             d_10 = interp1(pres,temp,10);
%             [a,b] = sort(temp);
%             mld(i,:) = interp1(a,pres(b),d_10-0.2);
%             clear mldidx; mldidx = find(flipud(pres) <= mld(i,:));
        elseif strcmp(method,'HT') == 1
            %Holte and Talley (density) method:
            mldindex = i; yesplot = 0; findmld
            mld(i,:) = mixeddp(mldindex); %From density
            clear mldidx; mldidx = find(flipud(pres) <= mld(i,:));
        else
            disp('Pick a method fool!')
        end
    end
    
    % Calculate mixed layer averages?
    if mean_calcs == 0
        % Nope
        mldmean = [];
        mldstd  = [];
    else
        % Calculate ML means for this profile of float data, loop through variables
        % What fields are in float structure?
        g = fields(f);
        for j = 1:length(fieldstr)
            clear varschar; varschar = char(fieldstr(:,j)); %find current variable name
            clear a;a = strmatch(fieldstr(j),g,'exact');
            if isempty(a) == 1
            else
                clear currvar; currvar = f.(varschar);
                if nansum(currvar(currentidx))==0 %address empty profiles
                    % Nope
                    mldmean.(varschar)(i,:) = 0;
                    mldstd.(varschar)(i,:)  = 0;
                elseif ~isempty(mldidx)
                    clear bad; bad = isnan(currvar(currentidx(mldidx)))==1;
                    currvar(currentidx(mldidx(bad))) = NaN;
                    if contains(varschar,'poc')
                        % If the variable is POC or bbp - use the median instead (Andrea suggested)
                        mldmean.(varschar)(i,:) = nanmedian(currvar(currentidx(mldidx)));
                    elseif contains(varschar,'bp700')
                        mldmean.(varschar)(i,:) = nanmedian(currvar(currentidx(mldidx)));
                    else % If the variable is not POC - take the mean
                        mldmean.(varschar)(i,:) = nanmean(currvar(currentidx(mldidx)));
                    end
                    mldstd.(varschar)(i,:) = nanstd(currvar(currentidx(mldidx)));
                else
                    % Set to surface value
                    clear tmp; tmp = find(~isnan(currvar(currentidx)));
                    clear surfaceidx; surfaceidx = tmp(end);
                    mldmean.(varschar)(i,:) = currvar(currentidx(surfaceidx));
                    mldstd.(varschar)(i,:)  = 0;
                end
            end
        end
        
        % Save date, profile, lat, and lon seperately to avoid being NaN in some cases
        for k = 1:length(fieldstr2)
            cc = char(fieldstr2(k)); %Current variable
            mldmean.(cc)(i,:) = nanmean(f.(cc)(currentidx)); %It's ok if it's the full profile mean for date
        end
        % Add MLD to the MLD mean float structure
        mldmean.mld = mld;
        % Save units for averaged variables
        n = 1;
        for j = 1:length(fieldstr)
            clear a;a = find(contains(g,fieldstr(j))==1, 1);
            if isempty(a) == 0
                desc.units(n,1) = fieldstr(j);
                desc.units(n,2) = units(j);
                
                desc.units(n,1) = fieldstr(j);
                desc.units(n,2) = units(j);
                n = n + 1;
            else
            end
        end
    end
end


