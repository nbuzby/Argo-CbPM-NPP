function extract_sat_data(bpath, fpath, satpath, movtspath, float_ids)
% extract_sat_data(bpath, fpath, satpath)
% finds and __ satellite data associated with float paths
%
%
% INPUT VARIABLES:
%   bpath:      base path
%   fpath:      directory of float files
%   satpath:    where composite satellite data is stored
%   movtspath:  where will save satellite data extracted at float location
%
% OPTIONAL INPUTS:
%   float_ids:  array of floats of interest, if not included will compute
%               for all existing downloaded float data
%
% OUTPUT 
%   structure 'uf' labeled by float ID of matching satellite data and
%   associated time/distance from each float profile
% 
% Functions called:
%   pos2dist

%%
var_str = {'PAR_par','CHL_chlor_a','GSM_bbp_443_gsm','KD490_Kd_490'}; %for file names
variables = {'par', 'chlor_a','bbp_443_gsm','Kd_490'}; %for folder paths

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

%establish limits to shorten ncread times
nclim = {'dim','min','max','length';'lon',2000,4500,4500-2000;'lat',2000,3500,3500-2000};

%% Extract & match sat data
for i=1:length(float_ids)
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = char(float_ids(i));
    else
        floatID = num2str(float_ids(i));
    end  
    load([fpath floatID '.mat'],'f') 

    uf.profile = f.profile(1,:);
    [prof,~] = unique(f.profile);
   
	for v=1:length(variables)
        cd(satpath)
        sat_file = ['Aqua.L3m_8D_' var_str{v} '_4km.nc'];
        disp(['Finding float ',floatID,' ',variables{v},' values'])
        
        %establish lat/lon bounding box for float to speed variable ncread
        lat = ncread(sat_file,'lat'); lat = lat([nclim{3,2}:nclim{3,3}-1]);
        lon = ncread(sat_file,'lon'); lon = lon([nclim{2,2}:nclim{2,3}-1]); %netcdf file is in deg E/W (-/+),
        time = ncread(sat_file, 'time')+7.1953e+05; %netcdf is in days since 1970-01-001
        
        for p = 1:length(prof)
            %find date of closest time-aligned composite file
            [duration,t_idx] = min(abs((time+4) - f.date(1,p))); %add 4 to time-center 8-day composite
            %save time difference
            uf.match.time.(variables{v})(p) = duration;
            
            uf.date(p) = time(t_idx);
            %larger matrix of satellite data from time-asligned composite
            var=flipud(ncread(sat_file,variables{v},[nclim{2,2},nclim{3,2},t_idx],[length(lon),length(lat),1])');
            
            if nansum(var,'all')== 0 %no data found within region
                match(p) = NaN;
                dist.(variables{v})(p) = NaN;
                uf.lat(p) = NaN;
                uf.lon(p) = NaN;
            else
                %% geographic match-up
                [lat_grid, lon_grid] = ndgrid(lat,lon);
                nandex = find(isnan(var));
                lon_grid(nandex) = NaN; lat_grid(nandex) = NaN;

                ref = table(var(~isnan(var)),lon_grid(~isnan(lon_grid)),lat_grid(~isnan(lat_grid)),...
                    'VariableNames',{'var_val','lon','lat'}); %eliminates NaN entries
                ref.closest = abs(ref.lon-f.lon(1,p))+abs(ref.lat-f.lat(1,p));
                [~,match_idx] = min(ref.closest);
                match(p) = ref.var_val(match_idx);
                %save distance difference in km
                dist.(variables{v})(p) = pos2dist(f.lat(1,p),f.lon(1,p),ref.lat(match_idx),ref.lon(match_idx),1);

                uf.lat(p) = ref.lat(match_idx);
                uf.lon(p) = ref.lon(match_idx);
            end
        end
        uf.(variables{v}) = match;
    end
    uf.match.distance = dist;
    cd(fpath)
	save([movtspath floatID '.mat'],'uf')
end

