%Find the 1st optical depth estimate and Zeu and save structures
%Based on mixed layer chlorophyll relationship and KdPAR
%
%Updated 24 July 2020 to calculated using kdPAR rather than kd490 and using
%estimated chla (LP regression with polyfit method) to calcaulte kd
% J Long (MBARI)

function [odmean,odstd,zeumean,zeustd] = calc_od_float(f,mldmean)

%Predefine struct variables
odmean = struct;
odstd = struct;
zeumean = struct;
zeustd = struct;

%Find unique profiles
[profile, profidx] = unique(f.profile);
%Possible variables and units... used to make a units field (notes) later
possible_vars = {'lon','lat','press','temp','sal','sigt','depth','oxy','oxysat','no3','chla','est_chla','qc_chla','LPCTD_adj_chla','bp700','bp700_corr','poc','ph_insitu','ph_25','ta_liar','dic_liar','co2_liar','cdom','sat_par','dens'};
units = {'Deg E','Deg N','dbar','Deg C','PSU','kg m^-3','m','umol kg^-1','%','umol kg^-1','ug L^-1','ug L^-1','ug L^-1','ug L^-1','m^-1','m^-1','mmol m^-3','total scale','total scale','umol kg^-1','umol kg^-1','uatm','ppb','E m^-2 d^-1','rho'};

%calculate density
if ~isfield(f,'dens')
    f.dens =  sw_dens(f.sal,f.temp,f.press); %kg/m3
else
end

for i = 1:length(profile)
    
    %Find data idx of current profile
    currentidx = find(f.profile == profile(i));
    
    if nansum(f.temp(currentidx))==0 %There's no data for this profile, move on
        OD(i,:) = NaN;
        Zeu(i,:) = NaN;
        zeuidx = '';
        odidx = '';
    else
        %% Find 1st Optical Depth
        %Get the ML-average chl for this profile
        calChlMLD = mldmean.chla(i); 
        Kd490_M07 = 0.0166 + 0.0773.*(calChlMLD).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
        Kd490_M07(Kd490_M07<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold (1986?), L&W Table 3.16
        KdPAR_M07 = 0.0665 + 0.874.*(Kd490_M07) - 0.00121./(Kd490_M07);  %Morel et al. 2007 Eq. 9' (KdPAR over 2 optical depths)
        OD(i,:) = 1/Kd490_M07;
        Zeu(i,:) = 4.6/KdPAR_M07;

        zeuidx = find(f.press(currentidx) <= Zeu(i,:));
        odidx = find(f.press(currentidx) <= OD(i,:)); %Be sure to flip the pressure back before finding index to match the rest of the variables
    end

    %% Calculate OD means for this profile of float data, loop through variables
    fieldstr = string(fieldnames(f)); %Get string of all field names
    for j = 1:length(fieldnames(f))
        varschar = char(fieldstr(j,:)); %find current variable name
        currvar = f.(varschar);
        if ~isempty(odidx) %HT method found the mld = 0 for at least one profile (f7601, I think profile 114) so needed to add this if statement
            bad = find(currvar(currentidx(odidx)) < 0); %None of these params should be negative
            currvar(currentidx(odidx(bad))) = NaN;
            if varschar(1:3) == 'poc' %It's POC => use median instead (Andrea suggested)
                odmean.(varschar)(i,:) = nanmedian(currvar(currentidx(odidx)));
                zeumean.(varschar)(i,:) = nanmedian(currvar(currentidx(zeuidx)));
            else %It's not POC => take the mean
                odmean.(varschar)(i,:) = nanmean(currvar(currentidx(odidx)));
                zeumean.(varschar)(i,:) = nanmean(currvar(currentidx(zeuidx)));
            end
            odstd.(varschar)(i,:) = nanstd(currvar(currentidx(odidx)));
            zeustd.(varschar)(i,:) = nanstd(currvar(currentidx(zeuidx)));
        else
            odmean.(varschar)(i,:) = NaN;
            odstd.(varschar)(i,:) = NaN;
            zeumean.(varschar)(i,:) = NaN;
            zeustd.(varschar)(i,:) = NaN;
        end
        
    end
end
%For some reason the HT method gives a NaN for certain rows, which
%makes the datetime = NaN, which is really hard to work with later on, so
%I'm going to remove those before saving
% clear nan_idx; nan_idx = find(~isfinite(odmean.date));
% if ~isempty(nan_idx)
%     disp('Found NaN')
%     for v = 1:length(fieldnames(f))
%         varschar = char(fieldstr(v,:)); %find current variable name
%         odmean.(varschar)(nan_idx) = NaN; %%% Just changed this from =[] to =NaN, which makes this obsolete, but I decided that I need to at least keep a row that represents this profile
%     end
% else
% end
% OD(nan_idx) = NaN;

zeumean.zeu = Zeu;
odmean.od = OD;

end