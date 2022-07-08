%% NPP MASTER CODE
clc; clear all; close all
initialize_argo()
do_pause()

%% Set paths - BASED ON COMPUTER SETUP
user = 'nbuzby';
pp = ['/raid/' user '/MATLAB/HOT-NPP-example-main/data/']; %path to build off of
%base path - the directory you will navigate from
bpath = [pp 'input/'];
%float data path - the directory you will save original float data to
fpath = [bpath 'float_data/Argo/'];
%depth averaged data path
dpath = [bpath 'depth_avg_data/'];
% Path to satellite data extracted at float location
movtspath = [bpath 'input_for_sat_comparison/Floats/'];
% Path to save NPP float data structure
npath = [pp 'npp/'];

%% Identify Floats - To find applicable floats in given region
cd(fpath)
lonlim = [-80 0];
latlim = [0 65];

% find floats with CHlA & BBP
[NAtl_profiles,NAtl_floats] = select_profiles(lonlim,latlim,[],[],'sensor','CHLA',...
    'outside','space');
[NAtl_profiles2,NAtl_floats2] = select_profiles(lonlim,latlim,[],[],'sensor','BBP700',...
    'outside','space');
good_float_ids = NAtl_floats2 .* ismember(NAtl_floats2, NAtl_floats);
idx = find(good_float_ids > 0);
good_float_ids = good_float_ids(idx);

display the number of matching floats and profiles
disp(' ');
disp(['# of matching profiles: ' num2str(length(NAtl_profiles))]);
disp(['# of matching floats: ' num2str(length(good_float_ids))]);

show_trajectories(good_float_ids,'color','multiple');

%% Step 1 - Process Float Data
% will look for what files you have already downloaded & processed and
% only download new float information

float_ids = good_float_ids;
variables = {'PRES','TEMP','PSAL','CHLA','BBP700'};

for i=1:length(float_ids)
    float_file = [num2str(float_ids(i)), '.mat'];
    if exist(float_file,'file') == 0 
        %download Argo SPROF files ('f' structures), formats and saves float data as mat file
        compile_argo(float_ids(i), variables, bpath, fpath) 
        
    elseif exist([num2str(float_ids(i)),'_depavg_data.mat'],'file')==0
        %average float over different depth horizons ('mld','od','zeu' structures)
        extract_depth_horizons(fpath, dpath, float_ids(i))
    else
        disp([num2str(float_ids(i)) ' has already been downloaded and depth-averaged'])
    end
end

%% Plot float parameters
% creates a figure for each float (within float_ids) showing Chlorophyll,
% Backscatter, and Float trajectory. Helpful to check for spiky profiles

for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([fpath floatID])
    
    figure(i)
    tiledlayout(1,3)
    nexttile()
    plot(f.chla,f.press,'-k')
    set(gca, 'YDir','reverse')
    title([floatID ' Chl-A'])
    nexttile()
    plot(f.bbp700,f.press,'-k')
    set(gca, 'YDir','reverse')
    title([floatID ' BBP 700'])
    nexttile()
    geoscatter(f.lat(1,:),f.lon(1,:),'.r')
    geobasemap grayland
end

%% Step 2 - Calculate NPP - WHERE MAY NEED TO CHANGE THINGS

%calculate depth-resolved npp from float data
npp_depth(dpath, fpath, movtspath, npath, float_ids)

%currently the way PAR is address in the function: (lines 111-113)
%%%   REPEATED FOR EACH FLOAT (specified by float_ids)
%     load([movtspath floatID '.mat']) %using float-matched satellite PAR values
%     clear irr; irr = uf.par(p); % 'uf.par' is a 1-by-# of profiles vector
%     clear uf 

% you could just save the par data for each float to the 'f' structures see
% PAR_datadownload side function

%% Depth-Resolved NPP Plots

%create plots for each float
for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([npath floatID '_CbPM.mat'])
    load([npath floatID '_sat_surface_CbPM.mat'])
    
    figure(i)
    tiledlayout(1,2)
    %Set Colortable
    cc = [repmat(rgb('white'),1,1);cmocean('tempo')];
    colormap(cmocean('tempo')); hold on;
    
    %Graph contour plot
    nexttile() 
    [~,hC] = contourf(f.date,f.press(:,1),f.npp_vec,50,':');hold on 
    plot(f.date,f.inpp.mld,'k','LineWidth',0.5);
    
    %contour plot aesthetics
    ylim([1 200])
    set(gca, 'YDir','reverse')
    datetick('x','mmm yyyy','KeepLimits','keepticks')
    ylabel('Depth (m)')
    hcb = colorbar; 
    titleString = 'NPP (mg C m^{-3} d^{-1})';
    r = ylabel(hcb,titleString);
    title(['NPP (CbPM) from ' floatID])
    
    %integrated NPP timeseries
    nexttile()
    plot(f.inpp.date, f.inpp.intNPP,'.-','color',[0.100342311936581,0.372347793683504,0.420205254896834],'linewidth',2); hold on
    plot(snpp.date(1,:),snpp.inpp.intNPP,'.-','Color',[0.470438935134762,0.692968336746152,0.555905872665848],'linewidth',2);
    datetick('x','mmm yyyy','KeepLimits','keepticks')
    title('Integrated NPP (200m)')
    ylabel('NPP (mg C m^{-3} d^{-1})')
    legend('Float (Integrated)','Satellite (Surface)')
    
    set(findall(gcf,'-property','Fontname'),'Fontname','Trebuchet','fontsize',17)
end

%% Integrated NPP Plots

tiledlayout(4,4) %change based on how many floats you're dealing with

for i=1:length(float_ids)
    floatID = num2str(float_ids(i));
    load([npath floatID '_CbPM.mat'])
    load([npath floatID '_sat_surface_CbPM.mat'])
    
    nexttile()
    plot(f.inpp.date, f.inpp.intNPP,'.-','linewidth',2); hold on
    datetick('x','mmm yyyy','KeepLimits','keepticks')
    title(['Float ' floatID])
    ylabel('NPP (mg C m^{-3} d^{-1})')
end
set(findall(gcf,'-property','Fontname'),'Fontname','Trebuchet','fontsize',17)

