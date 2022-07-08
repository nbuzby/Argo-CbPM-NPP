%% Setup
%set your base path - the directory you will navigate from
bpath = '/raid/nbuzby/MATLAB/HOT-NPP-example-main/';

%create and set your data path - the directory you will save data to
dpath = [bpath '/data/input/par_files/'];

cd(dpath)

%% Download files from OSU Website

MODISpar_url  = 'http://orca.science.oregonstate.edu/data/1x2/monthly/par.modis.r2018/hdf/';
year = [2012:2021];
month = {'001','032','060','061','091','092','121','122','152','153','182','183','213','214',...
    '244','245','274','275','305','306','335','336'};

for i =1:length(year)
    y = num2str(year(i));
    for z = 1:length(month)
        m = char(month(z));
        MODISpar_file = ['par.',y,m,'.hdf.gz'];

        from_str = [MODISpar_url,MODISpar_file]; % build target string
        to_str   = [dpath,MODISpar_file];    % build destination string

        [~,url_chk] = urlread(from_str); % See if file exisit on the web
        if url_chk == 1
            f = urlwrite(from_str,to_str);
            disp(' ');
            disp(['Data for float ',MODISpar_file,' retrieved from:  ',...
            MODISpar_url]);
            disp(['Saved as  ',f]);
        else
            disp(['No file found! For ',MODISpar_file])
        end
    end
end

%% Unzip files, delete, and get dates

cd(dpath)
gunzip('*.gz')
disp('Go delete the gz files!')

plist = dir('*.hdf');
par_files = {plist.name};
par_dates = struct('raw',strings(1,length(par_files)),...
    'serial_dateno',[]);

for i =1:length(par_files)
    par_file_char = char(par_files(i));
    par_dates.raw(1,i) = par_file_char(5:11);
    
    year = str2num(par_file_char(5:8));
    doy = str2num(par_file_char(9:11));
    par_dates.serial_dateno(i)= datenum(datetime(year,1,doy));
end

%% Read in and add PAR data (by location)
% ADDED TO 2A_get_npp_from_argo on 9/2/21
cd(dpath)

grid_lat = flip([-90:0.1667:90]); %par data has flipped layout
grid_lon = [-180:0.1667:180];

float_files = {'5904124.mat','5904172.mat'};
for i = 1:length(float_files)
    fl_char = char(float_files(i));
    floatID = fl_char(1:7);
    load([fpath fl_char])
    fl = f; clear f
    fl.par = [];
    
    [prof pidx] = unique(fl.profile);
    for p = 1:length(prof)
        current_idx = find(fl.profile == prof(p));
        
        %find closest date bewteen profile and PAR data
        [c idx_day] = min(abs(par_dates.serial_dateno-fl.date(p)));

        %read in associated HDF file (named by dates)
        n = find(contains(par_files, par_dates.raw(idx_day))==1);        
        filename = char(par_files(n));
        par_file = hdfinfo(filename);
        sds_info=par_file.SDS;
        par_data = hdfread(sds_info);

        lat = fl.lat(p);
        lon = 180-fl.lon(p); %lon values are provided in ºEast
        
        %find closest lat/lon values within the par data grid
        [latDistance, idx_lat] = min(abs(grid_lat-lat));
        [lonDistance, idx_lon] = min(abs(grid_lon-lon));
        
        %bring associated PAR value into float data structure
        fl.par(current_idx) = par_data(idx_lat, idx_lon);
        fl.par = fl.par.'; %transpose to column to match other parameters
    end
end

%% check float dates
% Pre 2021 Floats: 5904124,5904172

for i = 1:length(float_files)
    od_files_char = char(float_files(i));
    floatID = od_files_char(1:7);
    n = find(contains(f_files,floatID) == 1);
    fl_char = char(f_files(n));
    load([fpath fl_char]) %Load full profile data, f struct
  
    d_max = datestr(max(f.date));
    d_min = datestr(min(f.date));
    disp([floatID ' ' d_min ' to ' d_max])
    
end

