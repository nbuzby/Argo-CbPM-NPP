function compile_argo(float_ids, variables, bpath, fpath, plot_vars)
% compile_argo(float_ids, dpath, bpath)
% Specify floats and download/format associate data
%
% INPUT VARIABLES:
%   float_ids: float IDs of floats of interest
%   variables: parameters of interest
%   bpath: base path - the directory you will navigate from
%   fpath: data path - the directory you will save data to
%
% OPTIONAL INPUTS:
%   plot_vars: If want diagnostic plots of variables, provide cell array
%   (default will not make plots)
%
% OUTPUT 
%   will save .mat files of each associated float to specified dpath

%% translation keys 
%Key of variable names in ARGO and the associated output name herein
%These are the variables needed for NPP calcs
key1 = {'','TIME','LONGITUDE','LATITUDE','PRES','TEMP','PSAL','CHLA','BBP700','NITRATE',...
    'PARAMETER_DATA_MODE'}; 
%excluding cycle number as profile b/c sometimes have duplicate cycles when different lat/lon

%These are the field names you're going to save them as
key2 = {'','date','lon','lat','press','temp','sal','chla','bbp700','no3'};

%% Download Argo data

if isa(float_ids, 'cell')
    for i=1:length(float_ids)
        float_ids2(i) = str2num(float_ids{i});
    end
    clear float_ids; float_ids = float_ids2;
    clear float_ids2
end

good_float_ids = download_multi_floats(float_ids);

[Data, Mdata] = load_float_data(good_float_ids, variables,[],'interp_lonlat','no'); %exclude under ice (interpolated locations)
floats = fieldnames(Data);

%% Format data and save as structure

%filter out QC flagged values
[Data_full] = qc_filter(Data,{'PRES';'TEMP';'PSAL';'NITRATE'}, [1,2],...
    'CHLA',[1,2,5],'BBP700',[1,2,3],'POSITION',1);

for i=1:length(floats)
    fn = floats{i};
    floatID = fn(2:end);
    Data_good = Data_full.(floats{i}); %pair down to single float
    if min(min(Data_good.PRES)) < 0  
        xq = (0:1:1000)';
    else
        xq = (floor(min(min(Data_good.PRES))):1:1000)';%minimum depth of float to 1000m depth (to include max MLDs)
    end
    f = struct(); %empty structure for interpolated data
    
    %Format float data into structure with matched up variable names
    vars = fieldnames(rmfield(Data_good,'JULD'));  %don't need julian days field
    for j=1:length(vars)
        v_check = ismember(key1,vars(j));
        if all(v_check==0) %meaning none of the key terms matches the variable
        else
            v = find(ismember(key1,vars(j))==1);
            s.(key2{v}) = Data_good.(vars{j});
            f.(key2{v}) = [];
        end
    end
    
    %generate profile variable,can't always trust cycle number 
    rep_leng = length(xq);
    for n=1:(width(s.date))
        if n < width(s.date) %repeat until final entry
            if isequal(s.date(1,n),s.date(1,n+1)) && s.lat(1,n)== s.lat(1,n+1) && s.lon(1,n)== s.lon(1,n+1)
                disp([floatID ' duplicate profile ',num2str(x),' on ', datestr(s.date(1,n))]);
            else
               f.profile(:,n) = repelem(n,rep_leng); 
               s.profile(:,n) = repelem(n,height(s.date));
            end
        else %check backwards on final entry
            if isequal(s.date(1,n),s.date(1,n-1)) && s.lat(1,n)== s.lat(1,n-1) && s.lon(1,n)== s.lon(1,n-1)
                disp([floatID ' duplicate profile ',num2str(x),' on ', datestr(s.date(1,n))]);
            else
               f.profile(:,n) = repelem(n,rep_leng);
               s.profile(:,n) = repelem(n,height(s.date));
            end
        end
    end
    
    for v=2:length(vars)
        v_check = ismember(key1,vars(j));
        if all(v_check==0) %meaning none of the key terms matches the variable
        else
            %interpolate data onto 1m grid, for each profile
            prof = unique(f.profile);
            for p=1:length(prof)
                if any(strcmp({'temp','sal','chla','bbp700','no3'}, key2{v})) %only interpolate pressure-dependent data
                    % Get index of current profile
                    pidx = find(s.profile == prof(p));
                    % find fininte pressure values
                    idx = find(~isnan(s.press(pidx))); 
                    if isempty(idx)
                        f.(key2{v}) = [f.(key2{v}),NaN(length(xq),1)];
                    else
                        x = s.press(pidx); x= x(idx); %finite pressure values of given profile
                        val=s.(key2{v})(pidx); val=val(idx);

                        nandex = ~isnan(val); %find finite variable values
                        if sum(nandex) < 2 %empty profile, b/c need at least two points to interpolate
                            f.(key2{v}) = [f.(key2{v}),NaN(length(xq),1)];
                        else
                            intp_val = interpn(x(nandex),val(nandex),xq,'pchip');  %xq, float surface:1000m
                            %plot(val,-x,'o',intp_val,-xq,'-') %CHECK INTERPOLATION
                            f.(key2{v}) = [f.(key2{v}),intp_val];
                        end
                    end
                elseif strcmp(key2{v},'press')
                    f.press = repelem(xq,[1],[length(prof)]);
                else
                    %repeat profile data over new depth grid
                    f.(key2{v}) = repelem(s.(key2{v})(1,:),[length(xq)],[1]);
                end
            end
        end
    end
    save([fpath '/' floatID '.mat'],'-struct','f'); clear s f 
end

%% Plotting Checks
% assign default qc_flags if none provided as input
if nargin > 4
    for i=1:length(plot_vars)
        if strcmp(plot_vars{i},'CHLA')
            qc_flags = [1,2,5];
        else
            qc_flags = [1,2];
        end
        plot_profiles(Data, Mdata, plot_vars(i),'obs',1,'qc',qc_flags)
    end
end

