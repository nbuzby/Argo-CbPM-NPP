function extract_sat_data(bpath, fpath, satpath, movtspath, float_ids)
% extract_sat_data(bpath, fpath, satpath, movtspath)
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
%   thresholds: 2x1 array of grid point thresholds based on satellite 
%               autocorrelation analysis [vertical,horizontal] ~NOPE~ 
%
% OUTPUT 
%   structure 'uf' labeled by float ID of matching satellite data and
%   associated time/distance from each float profile
% 
% Functions called:
%   pos2dist

%%
var_str = {'chlor_a','par','GSM_bbp_443_gsm','Kd_490'}; %for file names
variables = {'chlor_a','par','bbp_443_gsm','Kd_490'}; %for folder paths
platform = {'uf','f'}; %PAR-based matchup should be saved to float structure, Chl-based matchup saved to sat structure
thresholds = [21,36];

if nargin >4 %are float ids or thresholds specified?
else
    %create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
    end 
end

lonlim = [-90 0]; 
latlim = [0 75];

%% Extract & match sat data
for i=1:length(float_ids)
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = char(float_ids(i));
    else
        floatID = num2str(float_ids(i));
    end  
    f=load([fpath floatID '.mat']); 

    uf.lat = f.lat(1,:);
    uf.lon = f.lon(1,:);
    uf.profile = f.profile(1,:);
    uf.date = f.date(1,:);
    [prof,~] = unique(f.profile);
    
    disp(['Finding float ',floatID,' matchups'])
   
    for v=1:length(variables)
        cd(satpath)
        if v==3 %bbp data have not undergone 2022 processing
            cd([satpath '/Release_2018.1QL'])
        end
        sat_file = ['Aqua.L3m_8D_' var_str{v} '_4km.nc'];
        
        %establish lat/lon bounding box for float to speed variable ncread
        lat = ncread(sat_file,'lat');
        [~,lat_idx1] = min(abs(lat-latlim(1))); 
        [~,lat_idx2] = min(abs(lat-latlim(2)));
        lat = flipud(lat(lat_idx1:lat_idx2-1)); %flip so array organized high-low 
        
        lon = ncread(sat_file,'lon'); 
        [~,lon_idx1] = min(abs(lon-lonlim(1))); 
        [~,lon_idx2] = min(abs(lon-lonlim(2)));
        lon = lon(lon_idx1:lon_idx2-1);

        [lat_grid, lon_grid] = ndgrid(lat,lon);
        deg_dist = abs(lat(2)-lat(1)); %distance in degrees between grid points

        time = ncread(sat_file, 'time')+719529; %netcdf is in days since 1970-01-001
        
        for p = 1:length(prof)
            if v==1 || v==2 %only matchup PAR & Chl
                %find date of closest time-aligned composite file
                [duration,t_idx] = min(abs((time+4) - f.date(1,p))); %add 4 to time-center 8-day composite
                %save time difference
                eval([platform{v} '.match.time(p) = duration;'])
                eval([platform{v} '.match.date(p) = time(t_idx);'])
    
                %larger matrix of satellite data from time-asligned composite
                var=flipud(ncread(sat_file,variables{v},[lon_idx1,lat_idx1,t_idx],[length(lon),length(lat),1])');
                
                if nansum(var,'all')== 0 || isnan(f.lat(1,p)) %no data found within region or profile not geolocated
                    eval([platform{v} '.' variables{v} '(p) = NaN;'])
                    eval([platform{v} '.match.distance.km(p) = NaN;'])
                    eval([platform{v} '.match.distance.vert(p) = NaN;'])
                    eval([platform{v} '.match.distance.horz(p) = NaN;'])
                    eval([platform{v} '.match.lat(p) = NaN;'])
                    eval([platform{v} '.match.lon(p) = NaN;'])
                else
                    %% geographic match-up
                    nandex = find(~isnan(var));
                    ref = table(var(nandex),lon_grid(nandex),lat_grid(nandex),...
                        'VariableNames',{'var_val','lon','lat'}); 
    
                    ref.closest = posdist(ref.lat,ref.lon,f.lat(1,p),f.lon(1,p),'s'); %great circle distance in km
                    ref.grid_lat = abs(ref.lat-f.lat(1,p))/deg_dist; %determine "grid-point" distance to work with autcorrelation threshold definition
                    ref.grid_lon = abs(ref.lon-f.lon(1,p))/deg_dist;

                    ref.closest(find(abs(ref.grid_lat) >thresholds(1))) = NaN; %vert autocorr threshold
                    ref.closest(find(abs(ref.grid_lon) >thresholds(2))) = NaN;%horz autocorr threshold

                    [val,match_idx] = min(ref.closest); 
                    if isnan(val) %no measurements within threshold
                        eval([platform{v} '.' variables{v} '(p) = NaN;'])
                        eval([platform{v} '.match.distance.km(p) = NaN;'])
                        eval([platform{v} '.match.distance.vert(p) = NaN;'])
                        eval([platform{v} '.match.distance.horz(p) = NaN;'])
                        eval([platform{v} '.match.lat(p) = NaN;'])
                        eval([platform{v} '.match.lon(p) = NaN;'])
                    else
                        eval([platform{v} '.' variables{v} '(p) = ref.var_val(match_idx);'])
                        
                        %save distance difference in km,grid points, and associated lat/lon positions
                        eval([platform{v} '.match.distance.km(p) = val;'])
                        eval([platform{v} '.match.distance.vert(p) = ref.grid_lat(match_idx);'])
                        eval([platform{v} '.match.distance.horz(p) = ref.grid_lon(match_idx);'])
                        
                        eval([platform{v} '.match.lat(p) = ref.lat(match_idx);'])
                        eval([platform{v} '.match.lon(p) = ref.lon(match_idx);'])
                    end
    
                end
            end
        end
        if v==3 || v==4 || v==2 %use Chl matchup for Kd,bbp,and sat-PAR
            ref_date = uf.match.date(p);
            t_idx = find(time==ref_date);

            var=flipud(ncread(sat_file,variables{v},[lon_idx1,lat_idx1,t_idx],[length(lon),length(lat),1])');
            %grid indices for matchup lat/lon
            [~,lat_idx] = min(abs(lat-uf.match.lat)); 
            [~,lon_idx] = min(abs(lon-uf.match.lon));
            %call associated variable values
            uf.(variables{v}) = var(sub2ind(size(var),lat_idx,lon_idx)); %sub2ind turns coordinate index values to linear index values
        end
    end
    disp('...done!')

	save([movtspath floatID '.mat'],'-struct','uf')
    save([fpath floatID '.mat'],'-struct','f','-append')
    clear f uf prof duration t_idx var match ref_lat ref_lon  nandex ref
end

