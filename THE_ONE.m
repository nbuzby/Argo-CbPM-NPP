%% NPP MASTER CODE
%Set paths
user = 'nbuzby';
ftype = 'Argo';
pp = ['/raid/' user '/MATLAB/HOT-NPP-example-main/data/']; %path to build off of
%base path - the directory you will navigate from
bpath = [pp 'input/'];
%float data path - the directory you will save original float data to
fpath = [bpath 'float_data/' ftype '/'];
%depth averaged data path
dpath = [bpath 'depth_avg_data/'];
%composite sat data path
satpath = '/raid/Data/Sat/MODIS/Aqua/Mapped/8-day/4km/'; 
% Path to satellite data extracted at float location
movtspath = [bpath 'input_for_sat_comparison/Floats/'];
% Path to save NPP float data structure
npath = [pp 'npp/'];

%variables = {'PRES','TEMP','PSAL','CHLA','BBP700','NITRATE'};
%cd(fpath); load('good_float_IDs.mat', 'float_ids');

%% Identify Floats
cd(fpath)

%establish NAtl polygon
lonlim = [-70,-80, -70, -15,  5, 30, -60]; 
latlim = [60, 30, 20,  20, 60, 80, 80];
%geoscatter(latlim,lonlim,30,'red','filled')

% find floats with CHlA & BBP
[NAtl_floats,NAtl_profiles] = select_profiles(lonlim,latlim,[],[],'sensor','DOWNWELLING_PAR',...
    'outside','space','interp_lonlat','no');
[NAtl_floats2,NAtl_profiles2] = select_profiles(lonlim,latlim,[],[],'sensor','BBP700',...
    'outside','space','interp_lonlat','no');
good_float_ids = NAtl_floats2 .* ismember(NAtl_floats2, NAtl_floats);
idx = find(good_float_ids > 0);
float_ids = good_float_ids(idx);

% display the number of matching floats and profiles
disp(' ');
%disp(['# of matching profiles: ' num2str(length(NAtl_profiles))]);
disp(['# of matching floats: ' num2str(length(float_ids))]);

show_trajectories(float_ids,'color','multiple');
clear NAtl_floats NAtl_floats2 NAtl_profiles NAtl_profiles2;

%% Step 1 - Process Float Data
variables = {'PRES','TEMP','PSAL','CHLA','BBP700','NITRATE','DOWNWELLING_PAR'};

cd(fpath)
n=1; m=1;
for i=1:length(float_ids)
    % try
    %     compile_argo(float_ids(i),variables,fpath,'DENS') 
    % catch
    %     bad.comp(n) = float_ids(i); n=n+1;
    % end 
    try
        cd(dpath)
        %1B average float over different depth horizons ('mld','od','zeu' structures)
        extract_depth_horizons(fpath, dpath, float_ids(i))
    catch
        bad.dep(m) = float_ids(i); m=m+1;
    end
end
    % float_file = num2str(float_ids(i));
    % cd(fpath)
    % compiled = isfile([float_file,'.mat']);
    % cd(dpath)
    % averaged = isfile([float_file,'_depavg_data.mat']);
    % 
    % if compiled == 0 %file does not exist
    %     cd(fpath)
    %     %1A download Argo SPROF files ('f' structures), formats and saves float data as mat file
    %     compile_argo(float_ids(i),variables,fpath) 
    %         %If want diagnostic plots add {'BBP700','CHLA'} to compile_argo function
    % end    
    % if averaged == 0 %float data has not been depth-averaged
    %     cd(dpath)
    %     %1B average float over different depth horizons ('mld','od','zeu' structures)
    %     extract_depth_horizons(fpath, dpath, float_ids(i))
    % else
    %     disp([num2str(float_ids(i)) ' has already been downloaded and depth-averaged'])
    % end
%end

%% Step 1C - Compile Satellite Data
lonlim = [-80 30]; latlim = [20 80];
%check to see if float-matched files exist
%for i=1:length(float_ids)
    cd(movtspath)
    % float_file = [num2str(float_ids(i)), '.mat'];
    % if isfile(float_file)
    %     %File exists. Are there any new data?
    %     load(float_file)
    %     file = dir(float_file);
    %     disp(['Float ', float_file(1:7), ' sat matchup was last modified ',file.date])
    % else
        %File doesn't exist already, run code
for i=1:length(float_ids)
            floatID = num2str(float_ids(i));
            % cd(movtspath)
            % %if ~isfile([floatID '.mat'])
            % try
            %     extract_sat_data(pp, latlim,lonlim,fpath, satpath, movtspath, float_ids(i))
            % catch
            %     issues{i,1} = ['issue with float ' num2str(floatID) ' sat matchup'];
            % end
            
        try
        %% float chl correction to sat values 
        threshold = 4;
        %for i=1:length(float_ids)
        floatID = num2str(float_ids(i));
        
        %f = load([fpath floatID '.mat']);
        load([movtspath floatID '.mat'],'chlor_a'); %sat chl data
        load([dpath floatID '_depavg_data.mat'],'odmean'); %od-averaged chl data
        %odmean = rmfield(odmean,{'gain','gain_flags','chl_corr'});

        odmean.gain.mean= mean(chlor_a,'omitnan')'./odmean.chla; %taking MEAN of rectangular region
        odmean.gain.med= median(chlor_a,'omitnan')'./odmean.chla; %taking MEDIAN of rectangular region
        odmean.gain_flags.mean = odmean.gain.mean.*0;
        odmean.gain_flags.med = odmean.gain.med.*0;

        odmean.gain_flags.mean(odmean.gain.mean > (mean(odmean.gain.mean,'omitnan')+iqr(odmean.gain.mean)*threshold))=1; %flag MEAN outliers
        odmean.gain_flags.med(odmean.gain.med > (mean(odmean.gain.med,'omitnan')+iqr(odmean.gain.med)*threshold))=1; %flag MED outliers


        odmean.chl_corr.mean = odmean.chla.*odmean.gain.mean; %NaNs out non-matches
        odmean.chl_corr.med = odmean.chla.*odmean.gain.med;
        save([dpath floatID '_depavg_data.mat'],'odmean','-append');
        
        floatID = num2str(float_ids(i));
        load([dpath floatID '_depavg_data.mat'], 'odmean')
        load([fpath floatID '.mat'], 'chla')
        chla_corr = chla.* odmean.gain.med';
        save([fpath floatID '.mat'],'chla_corr','-append'); clear chla chla_corr odmean
        %end
        catch
            issues{i,1} = ['issue with float ' num2str(floatID) ' odmean'];
        end
end

%% Step 2 - Calculate NPP

 for i=1:length(float_ids)
     floatID = num2str(float_ids(i));
     try
        npp_depth(dpath, fpath, npath, float_ids(i),'Y','sat')
    catch
        disp(['issue with ' num2str(float_ids(i)) ' depth npp, sat par, corrected']);
    end
    %  try
    %     npp_depth(dpath, fpath, npath, float_ids(i),'Y','prod')
    % catch
    %     disp(['issue with ' num2str(float_ids(i)) ' depth npp, ML par, corrected']);
    %  end
 %end
% %% 
    cd(npath)
    % try
    %     %if isfile([floatID,'_CbPM.mat'])==0 %float npp has not been calculated
    %     %2A calculate depth-resolved npp from float data
    %     npp_depth(dpath, fpath, npath, float_ids(i),'Y','sat')
    %     %end
    % catch
    %     disp(['issue with ' num2str(float_ids(i)) ' depth npp, sat par']);
    % end
    % try
    %     npp_depth(dpath, fpath, npath, float_ids(i),'Y','prod')
    % catch
    %     disp(['issue with ' num2str(float_ids(i)) ' depth npp, par product']);
    % end
    % try
    %     npp_depth(dpath, fpath, npath, float_ids(i),'N','prod')
    % catch
    %     disp(['issue with ' num2str(float_ids(i)) ' depth npp, uncorr par product']);
    % end
    try
        %if isfile([floatID,'_CbPM_uncorrected.mat'])==0
        npp_depth(dpath, fpath, npath, float_ids(i),'N')
        %end
    catch
        disp(['issue with ' num2str(float_ids(i)) ' uncorr depth npp']);
    end
    try
        %if isfile([floatID,'_sat_surface_CbPM.mat'])==0
        %2B calculate satellite npp SHOULD GO BEFORE 2A!!!!!
        npp_surface(pp, dpath, fpath, movtspath, npath, 'sat', float_ids(i),'mean','Y')
        %end
    catch
        disp(['issue with ' num2str(float_ids(i)) ' surf sat npp']);
    end
    try
        %if isfile([floatID '_float_uncorr_surface_CbPM.mat'])==0
        npp_surface(pp, dpath, fpath, movtspath, npath, 'float', float_ids(i),'mean','N')
        %end
    catch
        disp(['issue with ' num2str(float_ids(i)) ' surf float uncorr npp']);
    end
    % else
    %     disp([floatID ' NPP has been calculated at depth and surface'])
    %end
end


%% Step 2b - Integrate NPP

methds = {'_CbPM','_sat_surface_CbPM','_float_surface_CbPM','_CbPM_uncorrected','_float_uncorr_surface_CbPM'}; %,'_CbPM_parproduct','_CbPM_uncorrected_parproduct'};
%methds = {'_CbPM_parproduct'};

for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    for m=1:length(methds)
        load([dpath floatID '_depavg_data.mat'],'odmean','zeumean');
        zeu = zeumean.zeu; od = odmean.od; clear odmean zeumean
        load([fpath floatID '.mat'],'mld')
        if m==1 || m==4 %|| m==5 || m==6
            load([npath floatID methds{m} '.mat'],'press','inpp','npp_vec');
        else
            load([npath floatID methds{m} '.mat'],'depth','inpp','npp_vec')
            press = depth;
        end
        for p=1:width(press)
            [~,mld_idx] = min(abs(press(:,p)-mld(p)));
            inpp_mld(p) = sum(npp_vec(1:mld_idx,p));

            [~,zeu_idx] = min(abs(press(:,p)-zeu(p)));
            inpp_zeu(p) = sum(npp_vec(1:zeu_idx,p));
            [~,zeuplus_idx] = min(abs(press(:,p)-zeu(p)*1.5));
            inpp_zeuplus(p) = sum(npp_vec(1:zeuplus_idx,p));

            [~,od_idx] = min(abs(press(:,p)-od(p)));
            inpp_od(p) = sum(npp_vec(1:od_idx,p));
            inpp_under_od(p) = sum(npp_vec(od_idx+1:end,p));
            clear mld_idx zeu_idx zeuplus_idx od_idx
        end
        load([fpath floatID '.mat'],'lat','lon')
        lat = lat(1,:); lon=lon(1,:);
        inpp.inpp_mld = inpp_mld'; inpp.inpp_zeu = inpp_zeu'; inpp.inpp_zeuplus = inpp_zeuplus';
        inpp.inpp_od = inpp_od'; inpp.inpp_under_od = inpp_under_od;

        save([npath floatID methds{m} '.mat'],'-append','inpp','mld','zeu','od');
        clear press date_vec inpp npp_vec inpp_mld inpp_zeu inpp_zeuplus inpp_od inpp_under_od depth date lat lon
    end
end
%% Establish & Assign Fay & Mckliney Bioregions
%read in bioregions
cd([pp 'input'])
lonlim = [-90 15]; 
latlim = [0 90];

lat = ncread('Time_Varying_Biomes.nc','lat'); lon = ncread('Time_Varying_Biomes.nc','lon');
[~,lon_idx] = min(abs(lon-lonlim(1))); [~,lat_idx] = min(abs(lat-latlim(1)));
lon_count = abs(lonlim(1))+abs(lonlim(2));lat_count = abs(latlim(1))+abs(latlim(2));
lat = lat(lat_idx:(lat_idx+lat_count-1)); lon = lon(lon_idx:(lon_idx+lon_count-1));

bio_regions = ncread('Time_Varying_Biomes.nc','MeanBiomes',[lat_idx lon_idx],[lat_count lon_count]);
save('fay_mckinley_regions.mat','bio_regions','lat','lon')
%NOTE THE BIOME MATRIX IS FLIPPED, BUT SO ARE LAT/LON SO IT IS OKAY

% P = pcolor(fmr.lon,fmr.lat,fmr.bio_regions); P.EdgeColor = 'none'; cbh = colorbar; 
% hold on; xtickformat('degrees'); ytickformat('degrees')
% cbh.TickLabels = {'PEQU-E (6)','SP STPS (7)','NA ICE (8)','NA SPSS (9)','NA STSS (10)','NA STPS (11)','AEQU (12)','SA STPS (13)'};
% title('Fay & McKinley, 2014 Bioregions')

%assign biomes to profiles
fmr.biomes = ncread('biomes_2010_2022.nc','Biomes');
fmr.lat = ncread('biomes_2010_2022.nc','lat');
fmr.lon = ncread('biomes_2010_2022.nc','lon');
fmr.year = ncread('biomes_2010_2022.nc','year');

for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([fpath floatID '.mat'],'date','lat','lon');
    for p=1:length(lat(1,:))
        if ~isnan(date(1,p)) && ~isnan(lat(1,p))
            yr_match = year(date(1,p));
            if yr_match==2023 || yr_match==2024
                yr_match=2022;
            end
            tidx = find(fmr.year==yr_match);
            barray = fmr.biomes(:,:,tidx); barray(barray < 8 | barray > 11) = NaN; %only want NAtl regions  
            [~,idx2] = min(abs(fmr.lat - lat(1,p)));
            [~,idx1] = min(abs(fmr.lon - lon(1,p)));
            biome(p) = barray(idx1,idx2); %y=rows, x=columns
            
            if isnan(biome(p)) %find closest match
                [lon_grid, lat_grid] = ndgrid(fmr.lon,fmr.lat);
                nandex = isnan(barray);
                lat_grid(nandex) = NaN; lon_grid(nandex) = NaN;
            
                [~,tmp] = min(posdist(lat_grid(:),lon_grid(:),lat(1,p),lon(1,p)','s'));
                [~,lat_idx] = ismember(lat_grid(tmp),fmr.lat); 
                [~,lon_idx] = ismember(lon_grid(tmp),fmr.lon); 
                idx = sub2ind(size(barray),lon_idx,lat_idx);
                biome(p) = barray(idx);
                clear nandex tmp lat_idx lon_idx idx
            end
        else
            biome(p) = NaN;
        end
    end
    save([fpath floatID '.mat'],'-append','biome')
    clear date lat lon biome idx1 idx2 tidx biome barray
end

%separate out norwegian seas
poly = [56,5;69,-35;85,-90;85,60; 75,60;56,5];
%geoplot(poly(:,1),poly(:,2))

xv = poly(:,2); %lon coords
yv = poly(:,1); %lat coords
xq = tarray.Lon;
yq = tarray.Lat;

clear in on; [in,on] = inpolygon(xq,yq,xv,yv);

%apply to float files
ne_floats = unique(tarray.WMOID(in));
for i=1:length(ne_floats)
    floatID = num2str(ne_floats(i));
    f_idx = find(tarray.WMOID==ne_floats(i));
    biome = tarray.biome(f_idx);
    save([fpath floatID '.mat'],'biome','-append');
end
%% Combine all floats

wmoids = []; profs = []; months = []; lat = []; lon = []; date = []; intNPP = []; %delta = []; ...
    biomes = []; OD = []; MLD = []; Zeu = []; sat_inpp = []; fe_inpp = []; fu_inpp = []; 
    gain = []; gain_flag = []; Chl = []; Chl_sat = []; %inpp_mld =[]; inpp_zeu = []; inpp_od =[]; inpp_under_od = []; inpp_zeu_plus = [];
    %fp_inpp = []; fpu_inpp = [];
    
for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([fpath floatID '.mat'],'biome','mld')
    load([dpath floatID '_depavg_data.mat']);
    load([npath floatID '_CbPM.mat'],'inpp');
    fu = load([npath floatID '_CbPM_uncorrected.mat'],'inpp');
    s = load([npath floatID '_sat_surface_CbPM.mat'],'inpp');
    f = load([npath floatID '_float_surface_CbPM.mat'],'inpp');
    %fp = load([npath floatID '_CbPM_parproduct.mat'],'inpp');
    %fpu = load([npath floatID '_CbPM_uncorrected_parproduct.mat'],'inpp'); 

    %delta = [delta; abs(real(inpp.intNPP) - real(closest.inpp.intNPP)')./real(inpp.intNPP)];  %sat delta  
    intNPP = [intNPP; inpp.intNPP(:)];
    % inpp_mld = [inpp_mld; inpp.inpp_mld(:)];
    % inpp_zeu = [inpp_zeu; inpp.inpp_zeu(:)];
    % inpp_zeu_plus = [inpp_zeu_plus; inpp.inpp_zeuplus(:)];
    % inpp_od = [inpp_od; inpp.inpp_od(:)];
    % inpp_under_od = [inpp_under_od; inpp.inpp_under_od(:)];
    sat_inpp = [sat_inpp; s.inpp.intNPP'];
    fe_inpp = [fe_inpp; f.inpp.intNPP'];
    fu_inpp = [fu_inpp; fu.inpp.intNPP];
    %fp_inpp = [fp_inpp; fp.inpp.intNPP];
    %fpu_inpp = [fpu_inpp; fpu.inpp.intNPP];
    date = [date; odmean.date(1,:)'];
    lat = [lat; odmean.lat(1,:)']; lon = [lon; odmean.lon(1,:)'];
    months = [months;month(odmean.date(1,:))'];
    wmoids = [wmoids;repmat(str2num(floatID),length(odmean.date(1,:)'),1)];
    profs = [profs;odmean.profile(1,:)'];
    gain = [gain; odmean.gain.mean];
    gain_flag = [gain_flag; odmean.gain_flags.mean];
    MLD = [MLD; mld(1,:)'];
    OD = [OD; odmean.od'];
    Zeu = [Zeu; zeumean.zeu'];
    Chl = [Chl; odmean.chla];
    Chl_sat = [Chl_sat; odmean.chl_corr.med];

    if height(biome) > 1
        biomes = [biomes; biome];
    else
        biomes = [biomes; biome(:)];
    end


    clear mld odmean zeumean inpp f s biome fp fe fu fpu
end

tarray = table(wmoids, profs, lat, lon, date, months,biomes,OD,MLD,Zeu,gain,gain_flag,...
    Chl,Chl_sat,intNPP,sat_inpp,fe_inpp,fu_inpp,...%inpp_mld,inpp_zeu,inpp_od,inpp_under_od,inpp_zeu_plus,
    'VariableNames',... 
    {'WMOID','Prof','Lat','Lon','Date','Month','biome','od','mld','zeu','gain','gain_flag',...
    'chl_float','chla_sat','float_inpp','sat_inpp','float_extrap_inpp','float_uncorr_inpp'});
    %'float_inpp_mld','float_inpp_zeu','float_inpp_od', ...
    %'float_inpp_under_od','float_inpp_zeu_plus'});
    %fp_inpp,fpu_inpp,  'float_parprod_inpp','float_uncorr_parprod_inpp',
idx = tarray.Lon>90;
tarray.Lon(idx) = tarray.Lon(idx)-360;
clear wmoids profs lat lon date months intNPP delt biomes Zeu MLD OD fe_inpp sat_inpp ... 
fu_inpp gain gain_flag Chl Chl_sat fp_inpp fpu_inpp inpp_mld inpp_od inpp_zeu inpp_under_od inpp_zeuplus

cd('/raid/nbuzby/MATLAB/HOT-NPP-example-main/3_MAKE_FIGURES')
save('floats_deets.mat','tarray','-append')


var = []; 
for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    %load([dpath floatID '_depavg_data.mat'],'odmean')
    %load([fpath floatID '.mat'],'par');
    %load([movtspath floatID '.mat'],'chlor_a')
    load([npath floatID '_float_uncorr_surface_CbPM.mat'],'inpp')

    %var = [var; inpp.intNPP(:)];
    idx = find(tarray.WMOID == float_ids(i));
    tarray.float_extrap_uncorr(idx) = real(inpp.intNPP(1:length(idx)));
end
%tarray.float_extrap_uncorr = real(var);

%% Array of iNPP over different horizons
vars = {'float','sat','float_extrap','float_uncorr'}; %,'float_parprod'};
depths = {'mld','od','zeu','zeuplus','under_od'};

darray = struct();
for i=1:length(float_ids)
    floatID = num2str(float_ids(i));

    load([npath floatID '_CbPM.mat'],'inpp');
    fu = load([npath floatID '_CbPM_uncorrected.mat'],'inpp');
    s = load([npath floatID '_sat_surface_CbPM.mat'],'inpp');
    fe = load([npath floatID '_float_surface_CbPM.mat'],'inpp');
    %fp = load([npath floatID '_CbPM_parproduct.mat'],'inpp');
    for d=1:length(depths)
        for v=1:length(vars)
            input = ['inpp_' depths{d}];
            if i==1
                if v==1
                    darray.(input).(vars{v}) = inpp.(['inpp_' depths{d}])(:); 
                elseif v==2
                    darray.(input).(vars{v}) = s.inpp.(['inpp_' depths{d}])(:); 
                elseif v==3
                    darray.(input).(vars{v}) = fe.inpp.(['inpp_' depths{d}])(:); 
                elseif v==4
                    darray.(input).(vars{v}) = fu.inpp.(['inpp_' depths{d}])(:); 
                % elseif v==5
                %     darray.(input).(vars{v}) = fp.inpp.(['inpp_' depths{d}])(:);
                end
            else
                if v==1
                    darray.(input).(vars{v}) = [darray.(input).(vars{v}); inpp.(['inpp_' depths{d}])(:)]; 
                else
                    if v==2
                        tmp = s;       
                    elseif v==3
                        tmp = fe;
                    elseif v==4
                        tmp = fu;
                    % elseif v==5
                    %     tmp = fp;
                    end
                    darray.(input).(vars{v}) = [darray.(input).(vars{v}); real(tmp.inpp.(['inpp_' depths{d}])(:))];
                end
            end
        end
    end
    clear inpp fu s fe fp
end

for v=1:length(vars)
    darray.inpp_200.(vars{v}) = tarray.([vars{v} '_inpp']);
end
darray.inpp_200 = struct2table(darray.inpp_200);
darray.inpp_200.date = weeknum(tarray.Date);
darray.inpp_200.biome = tarray.biome;
for d=1:length(depths)
    darray.(['inpp_' depths{d}]) = struct2table(darray.(['inpp_' depths{d}]));
    darray.(['inpp_' depths{d}]).date = weeknum(tarray.Date);
    darray.(['inpp_' depths{d}]).biome = tarray.biome;
end
%darray.biome = tarray.biome; darray.date = tarray.date;

% cd('/raid/nbuzby/MATLAB/HOT-NPP-example-main/3_MAKE_FIGURES')
% save('floats_deets.mat','darray','-append')
%% Array of Depth-npp for various methods
T = table();
for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([npath floatID '_CbPM.mat'],'press','date_vec','npp_vec','mld','zeu','od'); 
    s = load([npath floatID '_sat_surface_CbPM.mat'],'npp_vec'); 
    fe = load([npath floatID '_float_surface_CbPM.mat'],'npp_vec');
    fu = load([npath floatID '_CbPM_uncorrected.mat'],'npp_vec');
    %fp = load([npath floatID '_CbPM_parproduct.mat'],'npp_vec');
    load([fpath floatID '.mat'],'biome','profile');

    f.press = press(:);
    %f.npp_delta_float = npp_delta.float(:);
    %f.npp_delta_sat = npp_delta.sat(:);
    f.npp_float = npp_vec(:);
    f.npp_sat = s.npp_vec(:); clear s 
    f.npp_extrp = fe.npp_vec(:); clear fe 
    f.npp_float_uncorr = fu.npp_vec(:); clear fu
    %f.npp_float_parprod = fp.npp_vec(:); clear fp;

    f.month = month(repmat(date_vec,200,1)); f.month = f.month(:);
    f.biome = repmat(biome(1:length(date_vec)),200,1); f.biome = f.biome(:);
    f.profile = repmat(profile(1,1:length(date_vec)),200,1); f.profile = f.profile(:);
    f.mld = repmat(mld(1,:),200,1); f.mld = f.mld(:);
    f.zeu = repmat(zeu,200,1); f.zeu = f.zeu(:);
    f.od = repmat(od,200,1); f.od = f.od(:);
    f.WMOID = repmat({floatID},width(press)*200,1);
    %f.npp_delta = 

    T = cat(1,T,struct2table(f));
    clear f press npp_delta date_vec biome profile inpp npp_vec fe fu s %fp
end
% cd('/raid/nbuzby/MATLAB/HOT-NPP-example-main/3_MAKE_FIGURES')
% save('floats_deets.mat','-append','T')

biomes = [8:11,14];
b_names = {'NA_ICE','NA_SPSS','NA_STSS','NA_STPS','Norwegian_Seas'};
vars = {'npp_float','npp_sat','npp_extrp','npp_float_uncorr'}; %,'npp_float_parprod'};
for v=1:length(vars)
    for b=1:length(biomes)
        idx = find(T.biome==biomes(b));
        P.(vars{v}).(b_names{b}) = pivot(T(idx,:),Columns='month',Rows='press',Method='mean',...
            DataVariable=vars{v});
        P_std.(vars{v}).(b_names{b}) = pivot(T(idx,:),Columns='month',Rows='press',Method='std',...
            DataVariable=vars{v});
        clear idx
    end
end

cd('/raid/nbuzby/MATLAB/HOT-NPP-example-main/3_MAKE_FIGURES')
save('float_deets.mat','-append','P','P_std')

%% Float Meta Data

%save sensor & DAC information to float files
cd(fpath); initialize_argo();
mf = struct(); 
for i=1:length(float_ids)
    %download_meta_files(str2num(float_ids{i})) 
    senpath = [fpath 'Meta/' float_ids{i} '_meta.nc'];
    sensor_ref = table(ncread(senpath,'SENSOR')',ncread(senpath,'SENSOR_MODEL')',...
        ncread(senpath,'SENSOR_MAKER')','VariableNames',{'Sensor','Model','Maker'});

    sensor_ref = table2cell(sensor_ref);
    idx = find(contains(sensor_ref,{'CHLA','BBP700'}));
    mf.id(i) = {floatID};
    mf.sensor(:,i) = sensor_ref(idx,1);
    mf.model(:,i) = sensor_ref(idx,2);
    mf.maker(:,i) = sensor_ref(idx,3);
    
    mf.center(i) = {ncread(senpath,'DATA_CENTRE')'};
    mf.inst(i) = {ncread(senpath,'OPERATING_INSTITUTION')'};
    clear sensor_ref idx senpath
end
save([fpath 'Meta/chl_sensor_meta_data_.mat'],'-struct','mf')

%using Argo Reference Table 4
% mf.center(find(strcmp(mf.center,'AO'))) = {'AOML, USA'};
% mf.center(find(strcmp(mf.center,'BO'))) = {'BODC, UK'};
% mf.center(find(strcmp(mf.center,'IF'))) = {'Ifremer, France'};

%sensor_options = unique(mf.model);


