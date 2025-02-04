function extract_depth_horizons(fpath, dpath, float_ids)

% extract_depth_horizons(fpath, dpath, float_ids)
%
% averages data over multiple depth horizons and save each, does not use using est_chla or Mari's bbp
%
% Functions called:
% calc_od_float.m
% calc_mld_float.m
%
% INPUT VARIABLES:
%   fpath: directory of float files
%   dpath: where will save depth-averaged files
%
% OPTIONAL INPUTS:
%   float_ids:  array of floats of interest, if not included will compute
%               for all existing downloaded float data
%
% OUTPUT 
%   will save each float averaged over multiple depth horizons as .mat files
%   Structures: mldmean, mldstd, odmean, odstd, zeumean, zeustd, and desc (units & info)


%% Find float data for averaging - based on input or what is already downloaded
%are float ids specified?
if nargin > 2
else %if not create list of existing files in downloaded float data directory
    cd(fpath)
    flist = dir('*.mat');
    float_files = {flist.name};
    for i=1:length(float_files)
        od_files_char = char(float_files(i));
        float_ids{i} = od_files_char(1:7);
    end
end

%% Get MLD averaged data first
for i = 1:length(float_ids)
    if isa(float_ids, 'cell') %floats are listed in cell array if taken from directory
        floatID = char(float_ids(i));
    else
        floatID = num2str(float_ids(i));
    end  
    f = load([fpath floatID '.mat']);
    
%     %Roesler correction
%     %['Chl was mult by two to remove the Roesler correction']};
%     f.chla = f.chla.*2;
%     chla =f.chla;
%     save([fpath '/' floatID '.mat'],'chla','-append')
    
    %%Extract depth averaged data
    [mldmean, mldstd] = calc_mld_float(f,'BY','N'); %using Holte & Talley MLD method
    
    %%Calculate OD and Zeu
    if f.lon(1)<0
        f.lon = f.lon+360; %od calculations don't work with negative values
    end
    [odmean,odstd,zeumean,zeustd] = calc_od_float(f,mldmean);
    if isfile([dpath '/' floatID '_depavg_data.mat'])
        save([dpath '/' floatID '_depavg_data.mat'],'mldmean','mldstd','odmean','odstd','zeumean','zeustd','-append')
    else
        save([dpath '/' floatID '_depavg_data.mat'],'mldmean','mldstd','odmean','odstd','zeumean','zeustd')
    end
    clear f mldmean mldstd zeumean zeustd odmean odstd
end