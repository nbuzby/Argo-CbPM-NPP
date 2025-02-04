function compile_argo(float_ids, variables, fpath, mld_method, plot_vars, OMZ)
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
%   mld_method: whether want MLD based on density of temperature threshold
%   ('DENS','TEMP') based on de Boyer Montégut et al. (2004). Default is
%   density
%   plot_vars: If want diagnostic plots of variables, provide cell array
%   (default will not make plots)
%   OMZ: specify Y/N on whether to incorporate oligotrophic/OMZ regions
%   dark correction (default will not use this function)
%
% OUTPUT 
%   will save .mat files of each associated float to specified dpath

%% optional variables
if nargin > 4
    a = exist('plot_vars');
    b = exist('OMZ');
    c = exist('mld_method');
    if a == 0
        plot_vars = 'N';
    elseif b == 0
        OMZ = 'N';
    elseif c==0
        mld_method = 'DENS';
    end
else
   plot_vars = 'N'; 
   OMZ = 'N';
   mld_method = 'DENS';
end

%% translation keys 
%Key of variable names in ARGO and the associated output name herein
%These are the variables needed for NPP calcs
key1 = {'','TIME','LONGITUDE','LATITUDE','PRES','TEMP','PSAL','CHLA_FINAL','BBP700_FINAL','NITRATE','MLD','DOWNWELLING_PAR'}; 
%excluding cycle number as profile b/c sometimes have duplicate cycles when different lat/lon

%These are the field names you're going to save them as
key2 = {'','date','lon','lat','press','temp','sal','chla','bbp700','no3','mld','irr'};

%% Download Argo data
cd(fpath)
if isa(float_ids, 'cell')
    for i=1:length(float_ids)
        float_ids2(i) = str2double(float_ids{i});
    end
    clear float_ids; float_ids = float_ids2;
    clear float_ids2
end


[Data, Mdata] = load_float_data(float_ids, variables,[],'interp_lonlat','no'); %exclude under ice (interpolated locations)
floats = fields(Data);

%saving basic variable QC info in meta files
tmp = {'PRES';'TEMP';'PSAL';'NITRATE';'DOWNWELLING_PAR'};
for f=1:length(floats)
    for t=1:length(tmp)
        if strcmp(tmp{t},'DOWNWELLING_PAR')
            if ~isfield(Data.(floats{f}),tmp{t})
                warning('float %s does not contain variable %s', floats{f}, ...
                    tmp{t});
            else
                if isfield(Data.(floats{f}),[tmp{t}, '_ADJUSTED']) && ...
                        sum(isfinite(Data.(floats{f}).([tmp{t}, ...
                        '_ADJUSTED'])(:)))
                    Mdata.(floats{f}).([tmp{t}, '_ADJUSTED_QC']) = Data.(floats{f}).([tmp{t}, '_ADJUSTED_QC']);
                else
                    warning('float %s adjusted values for %s are not available',...
                        floats{f}, tmp{t})
                    Data.(floats{f}) = rmfield(Data.(floats{f}),{tmp{t},[tmp{t} '_QC'],[tmp{t} '_dPRES'],...
                        [tmp{t} '_ADJUSTED'],[tmp{t} '_ADJUSTED_QC'],[tmp{t} '_ADJUSTED_ERROR'],['PROFILE_' tmp{t} '_QC']});
                end
            end
        else
            if ~isfield(Data.(floats{f}),tmp{t})
                warning('float %s does not contain variable %s', floats{f}, ...
                    tmp{t});
            else
                if isfield(Data.(floats{f}),[tmp{t}, '_ADJUSTED']) && ...
                        sum(isfinite(Data.(floats{f}).([tmp{t}, ...
                        '_ADJUSTED'])(:)))
                    Mdata.(floats{f}).([tmp{t}, '_ADJUSTED_QC']) = Data.(floats{f}).([tmp{t}, '_ADJUSTED_QC']);
                    if t==1 %save pressure data
                        Mdata.(floats{f}).([tmp{t}, '_ADJUSTED']) = Data.(floats{f}).([tmp{t}, '_ADJUSTED']);
                    end
                else
                    warning('float %s adjusted values for %s are not available, using raw values instead',...
                        floats{f}, tmp{t})
                    Mdata.(floats{f}).([tmp{t}, '_QC']) = Data.(floats{f}).([tmp{t},'_QC']);
                    if t==1 %save pressure data
                        Mdata.(floats{f}).(tmp{t}) = Data.(floats{f}).(tmp{t});
                    end
                end
            end
        end
        
    end
end

%% Pull QC-ed Data

%MLD CALC & basic parameter QC
for f=1:length(floats)
    Data.(floats{f}) = calc_auxil(Data.(floats{f}),'calc_dens',1,'calc_mld_dens',1,'calc_mld_temp',1);
end
[Data_full] = qc_filter(Data,{'PRES';'TEMP';'PSAL';'NITRATE'}, [1,2],'DOWNWELLING_PAR',[1,2,3,5,8]);

%% Opitcal Parameters QC

%CHLA(dark offset, range test, neg spike,NPQ)
%BBP (range, neg spike, offset test)
%PAR (dark identification, fitting to correct cloud/waves)
for f=1:length(float_ids)
    floatID = num2str(float_ids(f));
    %pair down to single float
    Data_raw = Data_full.(['F' floatID]);
    Data_raw.CHLA = Data.(['F' floatID]).CHLA;
    Data_raw.BBP700 = Data.(['F' floatID]).BBP700;
    
    %calc MLD
    if strcmp(mld_method,'DENS')
        Data_raw.MLD = Data.(['F' floatID]).MLD_DENS;
    elseif strcmp(mld_method,'TEMP')
        Data_raw.MLD = Data.(['F' floatID]).MLD_TEMP;
    end

    %PAR QC steps
    % sp = SolarSpectrum(lon, lat, date);
    
    for i=1:width(Data_raw.CYCLE_NUMBER)        
        %Chla Dark offset OMZ/Oligotrophic
        [pres, idx] = sort(Data_raw.PRES(:,i));
        chl = Data_raw.CHLA(idx,i);
        [uniqueVals, ~, idx] = unique(pres);
        %in situation of repeated measurements at same depths
        if length(pres) ~= length(uniqueVals)
            for u=1:length(uniqueVals)
                match_idx = find(idx==u);
                if length(match_idx)>1
                    pres(match_idx(2:end)) = NaN;
                    chl(match_idx(2:end)) = NaN;
                    chl(match_idx(1)) = median(chl(match_idx),'omitnan');
                end
            end
        end
        Data_raw.CHLA(:,i) = chl; Data_raw.PRES(:,i) = pres;
        Data_raw.CHLA_DARKOZ(:,i) = Darkoz(pres, chl, 0.5); 

        %Chla Dark offset
        if i==1
            [Data_raw.CHLA_DARK(:,i), Data_raw.NEWDARK(:,i), Data_raw.CHLA_QC(:,i)] = Dark(Data_raw.PRES(:,i), Data_raw.CHLA(:,i), ...
                [], Data_raw.MLD(i),i);
        elseif i>1
            [Data_raw.CHLA_DARK(:,i), Data_raw.NEWDARK(:,i), Data_raw.CHLA_QC(:,i)] = Dark(Data_raw.PRES(:,i), Data_raw.CHLA(:,i), ...
                Data_raw.NEWDARK(i-1), Data_raw.MLD(i),i);
        end
    
        if strcmp(OMZ, 'N') %north atlantic floats don't need OMZ method
            dark = 'CHLA_DARK';
        else
            dark = 'CHLA_DARKOZ';
        end
    
        %Range QC Tests
        [Data_raw.CHLA_RANGE(:,i), Data_raw.CHLA_QC(:,i)] = Range_test_chla(Data_raw.(dark)(:,i));
        [Data_raw.BBP700_RANGE(:,i), Data_raw.BBP700_QC(:,i)] = Range_test_bbp(Data_raw.BBP700(:,i));
    
        %Neg Spike Removal (QC=4)
        [Data_raw.CHLA_NEG(:,i), param_qc] = Neg_spike(Data_raw.(dark)(:,i));
        Data_raw.CHLA_QC(param_qc~=1,i) = param_qc(param_qc~=1);
    
        [Data_raw.BBP700_NEG(:,i), param_qc] = Neg_spike(Data_raw.BBP700(:,i));
        Data_raw.BBP700_QC(param_qc~=1,i) = param_qc(param_qc~=1);
    
        %Chla NPQ Correction
        %CHLA_NEG includes dark offset & negative spike removal (NOT RANGE?)
        Data_raw.CHLA_NPQ(:,i) = NPQ_MLD(Data_raw.CHLA_NEG(:,i),Data_raw.PRES(:,i),Data_raw.MLD(i));
    
        %REMOVE ROESSLER CORRECTION?
    
        %Bbp Offset
        bbp_qc = bbp_offset(Data_raw.BBP700(:,i));
        Data_raw.BBP700_QC(bbp_qc~=1,i) = bbp_qc(bbp_qc~=1);
    end
    Data_full.(['F' floatID])= Data_raw;
end
%% NEWER BBP MCQC FUNCTIONS
for f=1:length(float_ids)
    floatID = num2str(float_ids(f));
    Data_raw = Data_full.(['F' floatID]);
    if isfile([fpath 'Meta/' floatID '_meta.nc'])
        %Parking Depth meta data
        n = ncread([fpath 'Profiles/' floatID '_Sprof.nc'],'CONFIG_MISSION_NUMBER');
        m = ncread([fpath 'Meta/' floatID '_meta.nc'],'CONFIG_PARAMETER_NAME')'; m_idx = strcmp(m,"CONFIG_ParkPressure_dbar");
        p = ncread([fpath 'Meta/' floatID '_meta.nc'],'CONFIG_PARAMETER_VALUE');
    else
        n=NaN; m=NaN; p=NaN; m_idx =NaN;
    end

    for b={'missing','high_deep_val','neg_val','noisy'}%,'hook'}
        %bbp_missing: identify missing values
        %bbp_hihg_deep_val: find high values at depth
        %bbp_neg_val: identify negative values
        %bbp_noisy: flag profiles that are affected by noisy data
        
        for i=1:unique(Data_raw.CYCLE_NUMBER)
            if strcmp(b,'hook') %Save parking depth & check for "bad hook values"
                if length(p) < length(n)
                else
                    Data_raw.PARKING_DEPTH(i) = p(m_idx,n(i));
                    [bbp_qc] = bbp_hook(Data_raw.BBP700(:,i), Data_raw.PRES(:,i), Data_raw.PARKING_DEPTH(i));
                end
            else
                bbp_qc = eval(['bbp_' char(b) '(Data_raw.BBP700(:,i), Data_raw.PRES(:,i))']);
            end
            Data_raw.BBP700_QC(bbp_qc~=1,i) = bbp_qc(bbp_qc~=1);
        end
    end
    Data_full.(['F' floatID])= Data_raw;
end

%filter out QC flagged values
for f=1:length(float_ids)
    floatID = num2str(float_ids(f));
    sz = size(Data_full.(['F' floatID]).CHLA_QC);

    good_idx = find(Data_full.(['F' floatID]).CHLA_QC==1 | Data_full.(['F' floatID]).CHLA_QC==2 | Data_full.(['F' floatID]).CHLA_QC==5);
    Data_full.(['F' floatID]).CHLA_FINAL = NaN(sz);
    Data_full.(['F' floatID]).CHLA_FINAL(good_idx) = Data_full.(['F' floatID]).CHLA_NPQ(good_idx); clear good_idx;
    Data_full.(['F' floatID]).CHLA_FINAL(find(Data_full.(['F' floatID]).CHLA_FINAL < 0)) = NaN; %get rid of negatives

    good_idx = find(Data_full.(['F' floatID]).BBP700_QC==1 | Data_full.(['F' floatID]).BBP700_QC==2 | Data_full.(['F' floatID]).BBP700_QC==3);
    Data_full.(['F' floatID]).BBP700_FINAL = NaN(sz);
    Data_full.(['F' floatID]).BBP700_FINAL(good_idx) = Data_full.(['F' floatID]).BBP700_NEG(good_idx); clear good_idx;
    Data_full.(['F' floatID]).BBP700_FINAL(find(Data_full.(['F' floatID]).BBP700_FINAL < 0)) = NaN; %get rid of negatives
end

%save QC field info in metadata
for f=1:length(float_ids)
    floatID = num2str(float_ids(f));
    Mdata.(['F' floatID]).CHLA_QC = Data_full.(['F' floatID]).CHLA_QC;
    Mdata.(['F' floatID]).BBP700_QC = Data_full.(['F' floatID]).BBP700_QC;
end

%% Format QC-ed data and save as structure
clear f;
for i=1:length(floats)
    fn = floats{i};
    floatID = fn(2:end);
    Data_good = depth_interp(Data_full.(floats{i}),[],'prs_res',1);
    Data_good.MLD = Data_full.(floats{i}).MLD;
    
    %Format float data into structure with matched up variable names
    vars = fieldnames(Data_good); %don't need julian days field
    for j=1:length(vars)
        v_check = ismember(key1,vars(j));
        if all(v_check==0) %meaning none of the key terms matches the variable
        else
            v = find(ismember(key1,vars(j))==1);
            f.(key2{v}) = Data_good.(vars{j});
        end
    end
    
    %generate profile variable,can't always trust cycle number 
    for n=1:(width(f.date))
        if n < width(f.date) %repeat until final entry
            if isequal(f.date(1,n),f.date(1,n+1)) && f.lat(1,n)== f.lat(1,n+1) && f.lon(1,n)== f.lon(1,n+1) %same date & location
                disp([floatID ' duplicate profile ',num2str(n),' on ', datestr(f.date(1,n))]);
            else
               if n==1
                   f.profile = repelem(n,height(f.date))';
               else
                   f.profile(:,n) = repelem(n,height(f.date)); 
               end
            end
        else %check backwards on final entry
            if isequal(f.date(1,n),f.date(1,n-1)) && f.lat(1,n)== f.lat(1,n-1) && f.lon(1,n)== f.lon(1,n-1)
                disp([floatID ' duplicate profile ',num2str(n),' on ', datestr(f.date(1,n))]);
            else
               f.profile(:,n) = repelem(n,height(f.date));
            end
        end
    end

    if isfield(f,'irr')
        irr_int = NaN(size(f.irr));
        for p=1:width(f.irr)
            f.irr_int(:,p) = daily_PAR_given_depth(f.lon(1,p), f.lat(1,p), f.date(1,p),...
                f.irr(:,p), f.press(:,1));
        end
    end
    
    %save QC data to a metadata variable w/in structure
    nm = fields(Mdata.(floats{i})); idx = find(contains(nm,'QC') | contains(nm,'PRES'));
    for t=1:length(idx)
        f.meta.(nm{idx(t)}) = Mdata.(floats{i}).(nm{idx(t)});
    end
    f.meta.BBP700_QC = Data_full.(floats{i}).BBP700_QC;
    f.meta.CHLA_QC = Data_full.(floats{i}).CHLA_QC;
    
    cd(fpath)
    if isfile([fpath '/' floatID '.mat'])
        save([fpath '/' floatID '.mat'],'-struct','f','-append'); 
    else
        save([fpath '/' floatID '.mat'],'-struct','f');
    end
end
clear Data_full Data_good Data_raw 

 %% OLD INTERPOLATION CODE
% for i=1:length(floats)
%     fn = floats{i};
%     floatID = fn(2:end);
%     Data_good = Data_full.(floats{i}); %pair down to single float
% 
%     if min(min(Data_good.PRES)) < 0  
%         xq = (0:1:1000)';
%     else
%         xq = (floor(min(min(Data_good.PRES))):1:1000)';%minimum depth of float to 1000m depth (to include max MLDs)
%     end
%     f = struct(); %empty structure for interpolated data
% %     f.chla_mode_diff = sum(~isnan(Data.(floats{i}).CHLA_ADJUSTED)) - sum(~isnan(Data_good.CHLA));
% %     f.bbp_mode_diff = sum(~isnan(Data.(floats{i}).BBP700_ADJUSTED)) - sum(~isnan(Data_good.BBP700));
% 
%     %Format float data into structure with matched up variable names
%     vars = fieldnames(Data_good); %don't need julian days field
%     for j=1:length(vars)
%         v_check = ismember(key1,vars(j));
%         if all(v_check==0) %meaning none of the key terms matches the variable
%         else
%             v = find(ismember(key1,vars(j))==1);
%             s.(key2{v}) = Data_good.(vars{j});
%             f.(key2{v}) = [];
%         end
%     end
% 
%     %generate profile variable,can't always trust cycle number 
%     rep_leng = length(xq);
%     for n=1:(width(s.date))
%         if n < width(s.date) %repeat until final entry
%             if isequal(s.date(1,n),s.date(1,n+1)) && s.lat(1,n)== s.lat(1,n+1) && s.lon(1,n)== s.lon(1,n+1)
%                 disp([floatID ' duplicate profile ',num2str(n),' on ', datestr(s.date(1,n))]);
%             else
%                if n==1
%                    f.profile = repelem(n,rep_leng)';
%                    s.profile = repelem(n,height(s.date))';
%                else
%                    f.profile(:,n) = repelem(n,rep_leng); 
%                    s.profile(:,n) = repelem(n,height(s.date));
%                end
%             end
%         else %check backwards on final entry
%             if isequal(s.date(1,n),s.date(1,n-1)) && s.lat(1,n)== s.lat(1,n-1) && s.lon(1,n)== s.lon(1,n-1)
%                 disp([floatID ' duplicate profile ',num2str(n),' on ', datestr(s.date(1,n))]);
%             else
%                f.profile(:,n) = repelem(n,rep_leng);
%                s.profile(:,n) = repelem(n,height(s.date));
%             end
%         end
%     end
%     prof = f.profile(1,:);
%     for v=2:length(vars)
%         v_idx = ismember(key1,vars(v));
%         if all(v_idx==0) %meaning none of the key terms matches the variable
%         else
%             for p=1:length(prof)
%                 if any(strcmp({'temp','sal','chla','bbp700','no3'}, key2{v_idx})) 
%                     %only interpolate pressure-dependent data
%                     pidx = find(s.profile == prof(p)); % Get index of current profile
%                     % find fininte pressure values
%                     idx = find(~isnan(s.press(pidx))); 
%                     if isempty(idx)
%                         f.(key2{v_idx}) = [f.(key2{v_idx}),NaN(length(xq),1)];
%                     else
%                         x = s.press(pidx(idx)); %finite pressure values of given profile
%                         val=s.(key2{v_idx})(pidx(idx)); val(val==0)=NaN; 
%                         nandex = find(~isnan(val)); %find finite variable values
%                         if length(nandex) < 2 %empty profile, b/c need at least two points to interpolate
%                             f.(key2{v_idx}) = [f.(key2{v_idx}),NaN(length(xq),1)];
%                         else
%                             intp_val = interp1(x(nandex),val(nandex),xq,'pchip',NaN);  %xq, float surface:1000m
%                             %plot(val,-x,'o',intp_val,-xq,'-') %CHECK INTERPOLATION
%                             f.(key2{v_idx})(:,p) = intp_val;
%                         end
%                     end
%                 elseif strcmp(key2{v_idx},'press')
%                     f.press = repelem(xq,1,[length(prof)]);
%                 else
%                     %repeat profile data over new depth grid
%                     f.(key2{v_idx}) = repelem(s.(key2{v_idx})(1,:),[length(xq)],1);
%                 end
%             end
%         end
%     end
%     nm = fields(Mdata.(floats{i})); idx = find(contains(nm,'QC') | contains(nm,'PRES'));
%     for t=1:length(idx)
%         f.meta.(nm{idx(t)}) = Mdata.(floats{i}).(nm{idx(t)});
%     end
% 
%     if isfile([fpath '/' floatID '.mat'])
%         save([fpath '/' floatID '.mat'],'-struct','f','-append'); 
%     else
%         save([fpath '/' floatID '.mat'],'-struct','f');
%     end
%     clear s f 
% end
% clear Data_full Data_good Data_raw 

%% Plotting Checks
% assign default qc_flags if none provided as input
if nargin > 4
    if strcmp(plot_vars,'Y')
        for i=1:length(plot_vars)
            if strcmp(plot_vars{i},'CHLA')
                qc_flags = [1,2,5];
            else
                qc_flags = [1,2];
            end
            plot_profiles(Data, Mdata, plot_vars(i),'obs',1,'qc',qc_flags)
        end
    end
end

