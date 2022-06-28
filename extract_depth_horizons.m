function extract_depth_horizons(fpath, dpath, float_ids)

% extract_depth_horizons(fpath, dpath, float_ids)
%
% averages data over multiple depth horizons and save each, does not use using est_chla or Mari's bbp
%
% Functions called:
% find_numericID.m
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
    load([fpath floatID '.mat'],'f')
    
    %Roesler correction
    %['Chl was mult by two to remove the Roesler correction']};
    f.chla = f.chla.*2;
    
    %%Extract depth averaged data
    [mld, mldmean, mldstd, desc] = calc_mld_float(f,'HT',1,'Argo'); %using Holte & Talley MLD method
    nandex = find(isnan(mldmean.mld));
    %calc_mld fills NaN mld entries with 0 values --> replacing
    var_1=fieldnames(mldmean);
    for v=1:length(var_1)
        mldmean.(var_1{v})(nandex) = NaN;
    end
    var_2=fieldnames(mldstd);
    for v=1:length(var_2)
        mldstd.(var_2{v})(nandex) = NaN;
    end
    
    %%Calculate OD and Zeu
    if f.lon(1)<0
        f.lon = f.lon+360; %od calculations don't work with negative values
    end
    [odmean,odstd,zeumean,zeustd] = calc_od_float(f,mldmean);
    
    desc.info = {['MLD method used = HT'];...
        ['Kd used in Zeu and OD calculations = KdPAR'];...
        ['Date processed = ' char(datetime)];...
        ['poc_v2 calc using bbp_tot and Johnson et al 2017 relationship']};%...
    save([dpath '/' floatID '_depavg_data.mat'],'mldmean','mldstd','odmean','odstd','zeumean','zeustd','desc')
    
    clear f mldmean mldstd float_no zeumean zeustd vars odmean odstd
end