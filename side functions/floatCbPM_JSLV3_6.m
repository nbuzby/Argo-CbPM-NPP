%This function originally created by Meg Estapa (Skidmore), and accessed on
%August 27, 2019 when updates began by Jacki Long (MBARI), the original
%version by Estapa is "floatCbPM_MEstapa.m"
%
%V2 was created to be able to claculate NPP per profile rather than all at
%once, this modification was made because (unlike the float data in Estapa
%2019), the float data used herein may not by consistent in depth or depth
%resolution
%http://sites.science.oregonstate.edu/ocean.productivity/cbpm2.code.php
%
% Written by: Jacki Long (MBARI)
% JSL 08/28/2019 changed bbbasegrid443 to bbbasegrid700
% V2 Created 24 July 2020 (J Long) to use already-averaged OD data from
% "B_extract_data_at_depth_horizons.m"
%
% V3_4 was copied directly from "floatCbPM_JSL3_2_movTSdata.m" and adapted
% to use profile data of Chl and C as inputs
%
% V3_5 was copied from V3_4 to move one step further and define PAR, Ed,
% kd, with depth starting from the surface instead of using the median in
% the ML and then moving with depth after
%
% V3_6 was copied from V3_4 to remove any definitions based on being above
% or below the MLD, which is how the light field had been described
% previously. Now, PAR is attenuated starting from the surface and is
% progragated through the entire depth (similar to what Meg did)
%
% /* --------------------------- */
%
% THIS PROGRAM COMPUTES NPP W/ DEPTH RESOLVED DATA
%
% /* --------------------------- */

% /* --------------------------- */
% /*   input variables           */
% /* --------------------------- */
%       chl            chlorophyll concentration
%       bbp            backscatter
%       irr            Photosynthetically available radiation in Einsteins per
%                      day per square meter
%       k490           diffuse attenuation coefficient at 490 nm (units of 1 / m )
%       mld            mixing layer depth in meters
%       zno3           depth of the nitrocline
%       daylength      length of the day in decimal hours (calc from doy/posvec) 
%       *doy            day of year
%       *posvec         lattitude position

% /* --------------------------- */
% /*   depth resolved variables  */
% /* --------------------------- */
%
%    z[200];               /* depths */
%    chl_C[200];           /* chl:c ratio */
%    chlz[200];            /* chl */
%    mu[200];              /* growth */
%    Ezlambda[9][200];     /* fraction of light at nine wavelengths */
%    parz[200];            /* total light */
%    prcnt[200];           /* percent light */
%    Cz[200];              /* carbon */
%    ppz[200];             /* npp */
% /* --------------------------- */
% /*   output variables          */
% /* --------------------------- */  %
%    Primary productivity in milligrams Carbon per square meter per day
%    uMax;			/* max growth rate */
%    chlCarbonMax;		/* max chl:carbon ration */
%    nutTempFunc;		/* f(nut,T) */
%    chlCarbonSat;		/* satalite chl:carbon ratio */
%    carbon;		/* bbp converted to carbon */
%    IgFunc;		/* f(Ig) */
%    IgFuncz;               /* f(Ig) below the mixed layer depth */
%    z_eu;			/* euphotic depth at 1% light level */
%    npp;                   /* net primary production */

function [out] = floatCbPM_JSLV3_6(mld,Chl,bbp,k490,irr,posvec,doy,zvec,od)

daylength = day_length(doy,posvec(1));

load('/raid/nbuzby/MATLAB/HOT-NPP-example-main/2_CALC_NPP/functions/Mobley_Boss2012_PARconv.mat','PAR_air20min_mult','Ed_air20min_mult');

%% Here I'm basically just copying the code found at Oregon State's Productiviy Site for "Updated CbPM":
%    !Description:     opp_cbpm2 - computes daily primary productivity using a chl:Carbon ratio.
%                      This is a spectrally resolved version of the cbpm, using nine separate
%                      wavelengths.  It is also depth resolved, integrating the effects from
%                      the surface down to a fixed depth of 200 m.
%
%                      The cbpm2 algorithm estimates productivity using chl (m-1), bbp (m-1),
% 		     surface irradiance (Einsteins m-2 d-1), k490 (m-1), mld (m), zno3 (m)
%                      and day length (hours).
%
% Net primary productivity is carbon * growth rate, where carbon is proportional to particulate
% backscatter

% C = 13000 * (bbp - 0.00035);

% and growth rate is a function of nutrient and temperature stress (f(nut,T) and photoacclimation
% (f(Ig))
% 	mu = umax * f(nut,T) * f(Ig)
% where:
% 	umax = 2
% 	f(nut,T) = ((Chl/C)sat - y0) / ((Chl/C)max - y0)
% 	f(Ig) = 1 - exp (-5 * Ig)
% and:
% 	(Chl/C)sat = ratio of satellite observed chl and carbon (carbon from bbp)

% 	(Chl/C)max = 0.022 + (0.045-0.022) * exp (-3 * Ig)

% 	Ig = median mixed layer light level
% 	   = surface irradiance * exp (-k(lambda) * MLD/2)

%       The above items are analyzed for nine separate wavelengths, and is vertically resolved to a depth
% of 200 m.
%
% For more details, please see the paper by Westberry, et al (2008)
%
%    !Dependencies:
%       function austinPetzold_1986 ( double lambda, double K490 )
%
%          given a reference k490 vlaue, determine k(lambda) for a specified lambda
%
%          ref:
%             Austin, R. W., and T. J. Petzold (1986), Spectral dependence of the diffuse
%             attenuation coefficient of light in ocean waters, Opt. Eng., 25, 473 ? 479

%    !Revision History:

%    08-16-2010 first release version (Robert O'Malley)
%       [original code written in matlab by T. Westberry]
%
%    01-05-2011   O'Malley
%       add uMax trap on mu[m]
%       correct z_eu determination

%    !References and Credits
%
%       Westberry, T. Behrenfeld, M.J., Siegel, D.A., and Boss, E.; 2008.  Carbon-based
%       primary productivity modeling with vertically resolved photoacclimation.  Global
%       Biogeochemical Cycles, Vol. 22, GB2024, doi:10.1029/2007GB003078

% */

%    austinPetzold_1986( double, double );

% /* --------------------- */
% /*   spectral variables  */
% /* --------------------- */

lambda = [400, 412, 443, 490, 510, 555, 625, 670, 700];
parFraction = [0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];
Xcoeff = [.11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000];
eexp = [.64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000];
Kw= [.01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438];

%    Kd[9];
%    Kbio;
%    Kdif[9];
%
%    Klambda[9];
%    Eo[9];
%    Ez_mld[9];
%    par_mld;
%    delChlC;
%
%    y0;
%   int i;
%   int m;
%   int mzeu;
%    r;
%    prcnt0;
%    prcnt1;
%    z0;
%    z1;
%    numerator;
%    denominator;
%    fraction;
%   double deltaZ;
%
%   if(irr <= 0.0)
%     0.0
%   end

%   /* --------------------- */
%   /*   initialize values   */
%   /* --------------------- */

z_eu = -9999;    %//  1.05.2011
y0 = 0.0003;                     %/* min  chl:c  when  mu = 0 */
% z = [1:1:200]; %Have z be an array from 0 to 200 m
r = 0.1;

uMax = 2.0;                      %/* after BANSE (1991) */
npp = NaN; %Changed to NaN from 0.0
mzeu = 0;

%set up the wavelength-specific parts
for i=1:9
    Klambda(i) = AustinPetzold_CBpM2(lambda(i),k490); %This is meg's function
    Eo(i) = irr * parFraction(i);
    % Here we're saying delta z is mld/2 (median mld depth)
    Ez_mld(i) = Eo(i) * 0.975 * exp(-Klambda(i) * mld / 2.0);
end

%   /* ----------------------------- */
%   /*   reintegrate to get par at   */
%   /*   depth ...                   */
%   /*   do trapezoidal integration  */
%   /* ----------------------------- */

% Just realized I had i = 2:8, should have been 1:8
for i = 1:8
    par_mld(i) = (lambda(i+1)-lambda(i))*(Ez_mld(i+1)+Ez_mld(i))/2;
end
par_mld = sum(par_mld);

% Use this instead of loop above, results in ~0.02 higher value but I think
% it's more correct
% par_mld = trapz(lambda,squeeze(Ez_mld));% Integrating Edgrid with respect the to scalar spacing of lambda

% PAR is in E m-2 day-1 so need to divide by the fraction of actual daylight
% as a percentage to scale it
% This is used to calculate ChlC(max), so keeping it for now
par_mld = par_mld./daylength;

% Here we are using the median PAR in the MLD to define IgFunc (which if a
% variable that is later multiplied to get mu in the ML), so if I'm
% switching to depth-resolved Chl:C, it would make more sense probably to
% use PAR(z) instead of par_mld
% IgFunc = 1 - exp(-5.0 .* par_mld);

% ! I should actually change the way I use backscatter as well, here we are
% removing a "background of nonalgal particles", which is the 0.00035, but
% that should be similar to what we've removed in Mari's process... maybe?
% Her's is based on a background at depth, whereas this value is estimated
% from satellite obs (see Westberry et al 2008). But they are probably
% trying to get at this same background signal that we've already removed
% https://pscfiles.apl.uw.edu/woodgate/BeringStraitArchive/BeringStraitMooringData/BeringStraitMoorings2007to2009IPY_versionMar10/SeasoftForWavesProcessingSoftware/website/software/SBEDataProcforWindows.htm

% Find any low bbp values are change them
% Here is using bbp(440 or 443) relationship
bbp(bbp < 0.00035) = 0.00036;
carbon = 13000.0 .* (bbp - 0.00035);
% If using this relationship, the input bbp should be bbp(470) instead of
% 443
%carbon = 12128.*(bbp) + 0.59;  % scaling from Graff et al., 2015 only

% Find OD-averaged data
odidx = find(zvec < od);
ChlOD = nanmean(Chl(odidx));
carbonOD = nanmean(carbon(odidx));

chlCarbonSat = ChlOD / carbonOD;

if ( chlCarbonSat < y0 )
    chlCarbonSat = y0;
end

chlCarbonMax = 0.022 + (0.045-0.022) * exp(-3.0 * par_mld);
% % delChlC is used to calcaulte Chl:C later with depth, but not used in the
% % calculation of mu or NPP... so doesn't really matter
delChlC = chlCarbonMax - chlCarbonSat;

% This is used in the calculation of mu and might make sense to make it
% depth-dependent instead of based on chl:C satellite in the ML
% nutTempFunc = (chlCarbonSat - y0) ./ (chlCarbonMax - y0);

%   /* ''''''''''''''''''''''''' */
%   /*   calculate Kd offset     */
%   /*   carry through to depth  */
%   /*   non-chl attenuation     */
%   /* ------------------------- */

% Here we're calculating Kdif, which is used later in the Ed calc below the
% ML. At this point, I'm not sure about this, if it should stay as a
% relation to ChlOD or if I should modify it.
for i=1:9
    Kbio = Xcoeff(i) .* ChlOD .^eexp(i);
    Kd(i) = Kw(i) + Kbio;
    Kdif(i) = Klambda(i) - Kd(i);
end

%   /* ''''''''''''''''''''''''''''''''''' */
%   /*   integrate down the water column   */
%   /*   in one-meter steps                */
%   /* ----------------------------------- */

for m=1:200
    if m == 1
        % Describe the light field at just under the surface
        % This is how Meg did it, using these coefficients
        Edday = repmat([1:1:365],91,1);
        Edlat = repmat([90:-2:-90]',1,365); %Every 2 degrees
        %load MB12 coefficients to propagate Ed(0+) thru surface
        Edmult = interp2(Edday,Edlat,Ed_air20min_mult,doy,posvec(1));
        
        %         for i=1:9
        %             % This is how the original code propagated through surface, by
        %             % multiplying 0.975 and klamda, I geuss...
        Ezlambda(1,:) = Eo.*Edmult*exp(-Klambda(i)*zvec(m));
        %         end
    else
        for i=1:9
            Kbio = Xcoeff(i) * Chl(m-1) ^ eexp(i);     %/*  after Morel & Maritorena (2001)  */
            Kd(i) = Kw(i) + Kbio + Kdif(i);         
            clear del_z; del_z = abs(zvec(m) - zvec(m-1));
            Ezlambda(m,i) = Ezlambda(m-1,i)*exp(-Kd(i)*del_z);
        end
    end
    parz(m) = 0.0;
    for i = 1:8
        parz(m) = parz(m) + (lambda(i+1)-lambda(i))*(Ezlambda(m,i+1)+Ezlambda(m,i))/2;
    end
    
    IgFuncz = 1 - exp(-5.0 .* (parz(m)/daylength));
    
    chl_C(m) = Chl(m)/carbon(m);
    
    nutTempFunc = (chl_C(m) - y0) ./ (chlCarbonMax - y0);
    
    % Calculate growth rate
    mu(m) = uMax * nutTempFunc * IgFuncz;
    if mu(m) > uMax       
        mu(m) = uMax;       
    else; end            
        
    prcnt(m) = parz(m) / (irr * 0.975);
    
    %     /*  track this to get to the euphotic depth  */
    
    if prcnt(m) >= 0.01 %If the light available is > or equal to 1% then just record current m
        mzeu = m;
    else
        
        %       /* ''''''''''''''''''''''''''' */
        %       /*   now find 1% light depth   */
        %       /*   in case the user wants    */
        %       /*   to use this information   */
        %       /* --------------------------- */
        
        if z_eu == -9999     %Then it is the first time z_eu has been found, continue.  // 01.05.11
            prcnt0 = prcnt(mzeu); %Use the depth found on previous iteration
            prcnt1 = prcnt(mzeu+1); %Percent at one meter deeper
            z0 = zvec(mzeu);
            z1 = zvec(mzeu+1);
            numerator = prcnt0 - 0.01;
            denominator = prcnt0 - prcnt1;
            fraction = numerator / denominator;
            z_eu = z0 + (z1-z0)*fraction;
        end
        
    end
    ppz(m) = mu(m) * carbon(m);
    
end
%
%   /* ------------------------------- */
%   /*   do trapezoidal integration    */
%   /*   from m = 0 to m = 200         */
%   /* ------------------------------- */
%
%   //  note:  186 m is the euphotic depth for pure water

% Calcualte iNPP
% First set neg NPP to 0
ppz(ppz<0) = 0;
goodidx = find(~isnan(ppz));
out.inpp = trapz(zvec(goodidx),squeeze(ppz(goodidx)));% Integrating Edgrid with respect the to scalar spacing of lambda

%Transform some of the variables to save as an output array
out.ppz = ppz';
out.mu = mu';
out.chl_C = chl_C';
out.parz = parz';
out.z = zvec';
out.prcnt = prcnt';
out.carbon = carbon';
end




