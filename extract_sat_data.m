function extract_sat_data(pp,latlim,lonlim, fpath, satpath, movtspath, float_ids)
% extract_sat_data(bpath, fpath, satpath, movtspath)
% finds and pulls satellite data associated with float paths
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

%% Setup
cd([pp '/input/input_for_sat_comparison'])
zonal = ncread('sat_chl_variability.nc','zonal_med'); 
z_lat = ncread('sat_chl_variability.nc','z_lat'); z_lon = ncread('sat_chl_variability.nc','z_lon');
meridional = ncread('sat_chl_variability.nc','meridional_med'); 
m_lat = ncread('sat_chl_variability.nc','m_lat'); m_lon = ncread('sat_chl_variability.nc','m_lon');
variables = {'chlor_a','par','bbp_443','Kd_490'}; %for folder paths

if nargin >3 %are float ids or thresholds specified?
else
    %create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
        float_ids{i} = od_files_char(1:7);
    end 
end

%Establish *filled* interpolant products from surfaces for matchup distance lookup
[x,y] = ndgrid(m_lon,m_lat); Fy = griddedInterpolant(x,y,fillmissing2(meridional,'nearest'),'linear','none');
[x,y] = ndgrid(z_lon,z_lat); Fx = griddedInterpolant(x,y,fillmissing2(zonal,'nearest'),'linear','none');
clear x y

%% Extract & match sat data
%universal sat data dimensions
cd(satpath)
lat = ncread('Aqua.L3m_8D_chlor_a_4km.nc','lat'); 
lon = ncread('Aqua.L3m_8D_chlor_a_4km.nc','lon');
time = ncread('Aqua.L3m_8D_chlor_a_4km.nc', 'time')+719529; %netcdf is in days since 1970-01-001

% for v=1:length(variables)
%     time = ncread(['Aqua.L3m_8D_' variables{v} '_4km.nc'], 'time')+719529; %netcdf is in days since 1970-01-001
%     disp([variables{v} ' time length = ' num2str(length(time)) ' ' datestr(time(end))])
% end

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
    for v=1:length(variables)
        uf.(variables{v}) = NaN(size(uf.profile));
        if v==1
            uf.match.lat = NaN(size(uf.profile));
            uf.match.lon = NaN(size(uf.profile));
            uf.match.Ly = NaN(size(uf.profile));
            uf.match.Lx = NaN(size(uf.profile));
        end
    end
    
    disp(['Finding float ',floatID,' matchups'])
    %fprintf(1,['Finding float ',floatID,' matchups'],prog);
   
    for p = 1:length(prof)
        %find date of closest time-aligned composite file
        [duration,t_idx] = min(abs((time+4) - f.date(1,p))); %add 4 to time-center 8-day composite
        %save time difference
        uf.match.time(p) = duration;
        uf.match.date(p) = time(t_idx);

        if ~isnan(f.lat(1,p)) && ~isempty(t_idx) %only matchup for geolocated/time aligned profiles
            %% geographic match-up
            %get rectangle region sat_chl surfaces            
            Ly = floor(Fy(f.lon(1,p),f.lat(1,p))); %meridional surface lookup
            Lx = floor(Fx(f.lon(1,p),f.lat(1,p))); %zonal surface lookup
            uf.match.Ly(p) = Ly;
            uf.match.Lx(p) = Lx;

            if ~isnan(Ly) && ~isnan(Lx) %need finite distances
                %find associated indices
                [~,lat_idx] = min(abs(lat-f.lat(1,p)));
                [~,lon_idx] = min(abs(lon-f.lon(1,p)));
                [lat_grid,lon_grid] = ndgrid(lat(lat_idx-Ly:lat_idx+Ly),lon(lon_idx-Lx:lon_idx+Lx));
    
                for v=1:length(variables)
                    %read in sat data from region
                    var = ncread(['Aqua.L3m_8D_' variables{v} '_4km.nc'],variables{v},[lon_idx-Lx lat_idx-Ly t_idx],[2*Lx+1 2*Ly+1 1]); 
                    if sum(~isnan(var),'all') > 0 %have >0 measurement within threshold
                        uf.(variables{v})(1:((2*Ly+1)*(2*Lx+1)),p) = var(:);
                        if v==1 %only need to save locations once
                            %save associated lat/lon positions
                            uf.match.lat(1:((2*Ly+1)*(2*Lx+1)),p) = lat_grid(:);
                            uf.match.lon(1:((2*Ly+1)*(2*Lx+1)),p) = lon_grid(:);
                        end
                    end
                end
            end
        end
    end
    %nan out zeros from varying length arrary addtions
    for v=1:length(variables)
        uf.(variables{v})(uf.(variables{v})==0) = NaN;
    end
    uf.match.lat(uf.match.lat==0) = NaN;
    uf.match.lon(uf.match.lon==0) = NaN;
    if height(uf.par>1)
        f.par = mean(uf.par,'omitnan'); %save to regular float structure for npp calcs, and precalc mean
    end
    
    disp('...done!')
    save([movtspath floatID '.mat'],'-struct','uf')
    save([fpath floatID '.mat'],'-struct','f','-append')

    clear f uf prof duration t_idx var lat_idx lon_idx Ly Lx lat_grid lon_grid
end

