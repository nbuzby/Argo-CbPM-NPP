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
% Sept 12, 2020: just realized I wasn't setting neg NPP values to 0, but
% shouldn't have made too much of a diffference?
%
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
%       daylength      length of the day in decimal hours.
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


function [out] = floatCbPM_JSLV3_2_movTSdata(mld,Chl,bbp,k490,zno3,irr,daylength)

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
z = [1:1:200]; %Have z be an array from 0 to 200 m
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
% for i = 1:8
%     par_mld(i) = (lambda(i+1)-lambda(i))*(Ez_mld(i+1)+Ez_mld(i))/2;
% end
% par_mld = sum(par_mld);

% Use this instead of loop above, results in ~0.02 higher value but I think
% it's more correct
 par_mld = trapz(lambda,squeeze(Ez_mld));% Integrating Edgrid with respect the to scalar spacing of lambda

% PAR is in E m-2 day-1 so need to multiply by the fraction of actual daylight
% as a percentage to scale it
%   par_mld /= daylength;
% Here I had been dividing PAR by 24, but I think it's supposed to be the
% actual length of daylight hours (in "decimal hours")
par_mld = par_mld./daylength;

IgFunc = 1 - exp(-5.0 .* par_mld);

if(bbp < 0.00035)
    bbp = 0.00036;
end
carbon = 13000.0 * (bbp - 0.00035);

chlCarbonSat = Chl / carbon;

if ( chlCarbonSat < y0 )
    chlCarbonSat = y0;
end

chlCarbonMax = 0.022 + (0.045-0.022) * exp(-3.0 * par_mld);
delChlC = chlCarbonMax - chlCarbonSat;

nutTempFunc = (chlCarbonSat - y0) ./ (chlCarbonMax - y0);

%   /* ''''''''''''''''''''''''' */
%   /*   calculate Kd offset     */
%   /*   carry through to depth  */
%   /*   non-chl attenuation     */
%   /* ------------------------- */

for i=1:9
    Kbio = Xcoeff(i) .* Chl .^eexp(i);
    Kd(i) = Kw(i) + Kbio;
    Kdif(i) = Klambda(i) - Kd(i);
end

%   /* ''''''''''''''''''''''''''''''''''' */
%   /*   integrate down the water column   */
%   /*   in one-meter steps                */
%   /* ----------------------------------- */

for m=1:200
    
    %     /* ---------------------------------------------- */
    %     /*   if you are in the mixed layer, do this way   */
    %     /* ---------------------------------------------- */
    %% In the mixed layer
    if z(m) < mld 
        chl_C(m) = chlCarbonSat; %chlCarbonSat = Chl / carbon
        chlz(m) = chl_C(m) * carbon;
        %nutTempFunc is dependent on chlCarbonSat
        mu(m) = uMax * nutTempFunc * IgFunc;
        
        if mu(m) > uMax       %//  1.05.2011
            mu(m) = uMax;           %//  1.05.2011
        else; end             %//  1.05.2011
        
        for i=1:9
            Ezlambda(m,i) = Eo(i)*0.975*exp(-Klambda(i)*z(m));
        end
        
        % I think this is where you're defining the growth irradiance as
        % the median?
        parz(m) = 0.0;
%         for i = 1:8
%             parz(m) = parz(m) + (lambda(i+1)-lambda(i))*(Ezlambda(m,i+1)+Ezlambda(m,i))/2;
%         end
        % Use this instead of loop
         parz(m) = trapz(lambda,squeeze(Ezlambda(m,:)));% Integrating Edgrid with respect the to scalar spacing of lambda
        
        Cz(m) = carbon;
        
    else
        %% You are deeper than the mixed layer
        %
        %       /* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
        %       /*   if below mixed layer must treat properties differently   */
        %       /* ---------------------------------------------------------- */
        
        for i=1:9
            Kbio = Xcoeff(i) * chlz(m-1) ^ eexp(i);     %/*  after Morel & Maritorena (2001)  */
%             Kd[i] = Kw[i] + Kbio;
%             Kd[i] += Kdif[i];
            Kd(i) = Kw(i) + Kbio + Kdif(i); % The original code (above), I think should be like this
            Ezlambda(m,i) = Ezlambda(m-1,i)*exp(-Kd(i)*1.0);
        end
        
        parz(m) = 0.0;
        for i = 1:8
            parz(m) = parz(m) + (lambda(i+1)-lambda(i))*(Ezlambda(m,i+1)+Ezlambda(m,i))/2;
        end
        % Use this instead of loop above
%         parz(m) = trapz(lambda,squeeze(Ezlambda(m,:)));% Integrating Edgrid with respect the to scalar spacing of lambda
        
        
        %      !!!!! Need depth of the nitrocline; from
        %       Westberry "?Nitracline depths (zNO3) were calculated from monthly
        %       climatological nutrient fields reported in the World Ocean Atlas
        %       [Conkright et al., 2002] and defined as the depth where nitrate + nitrite exceeded 0.5 uM.
        deltaZ = zno3 - z(m);
        if deltaZ < 0
            deltaZ = 0;
        else; end
        
        chl_C(m) = 0.022 + (0.045-0.022) * exp(-3.0 .* (parz(m) / daylength));%Westberry et al. 2008; Max possible under nut-replete conditions,
        % In this HNLC region, deltaZ is almost always = 0
        chl_C(m) = chl_C(m) - delChlC .* (1-exp(-0.075*deltaZ));

        IgFuncz = 1 - exp(-5.0 .* (parz(m)/daylength));
        mu(m) = uMax * nutTempFunc * IgFuncz;
        
        if mu(m) > uMax     %//  1.05.2011
            mu(m) = uMax;           %//  1.05.2011
        else; end              %//  1.05.2011
        
        if (mu(m-1) >= r )
            Cz(m) = carbon;
        else
            Cz(m) = carbon * mu(m-1) / r;
        end
        
        chlz(m) = chl_C(m) * Cz(m);
        
    end
    
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
            z0 = z(mzeu);
            z1 = z(mzeu+1);
            numerator = prcnt0 - 0.01;
            denominator = prcnt0 - prcnt1;
            fraction = numerator / denominator;
            z_eu = z0 + (z1-z0)*fraction;
        end
        
    end
    ppz(m) = mu(m) * Cz(m);
    
end
%
%   /* ------------------------------- */
%   /*   do trapezoidal integration    */
%   /*   from m = 0 to m = 200         */
%   /* ------------------------------- */
%
%   //  note:  186 m is the euphotic depth for pure water

% if mzeu < 186 %saying "if mzeu is greater than the euphotic depth of pure water"
%     %Calcualte iNPP
%     npp(1) = 0;
%     for i = 2:199
%         out.inpp = (z(i+1)-z(i))*(ppz(i+1)+ppz(i))/2;
%     end
% else
%     npp = -9999;
% end

% Use this instead of loop above
ppz(ppz<0) = 0;
out.inpp = trapz(z,squeeze(ppz));% Integrating Edgrid with respect the to scalar spacing of lambda

%Transform some of the variables to save as an output array
out.ppz = ppz';
out.mu = mu';
out.chlz = chlz';
out.Cz = Cz';
out.chl_C = chl_C';
out.parz = parz';
out.z = z';
out.prcnt = prcnt';
out.carbon = carbon';
end


%% /* =================================================================  */
%
% austinPetzold_1986 ( double lambda,
% double K490 )
%
% wave() =  350, 360, 370, 380, 390, 400,
% 410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
% 510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
% 610, 620, 630, 640, 650, 660, 670, 680, 690, 700 ;
%
% M() =  2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383,
% 1.7591, 1.6974, 1.6108, 1.5169, 1.4158, 1.3077, 1.1982, 1.0955, 1.0000, 0.9118,
% 0.8310, 0.7578, 0.6924, 0.6350, 0.5860, 0.5457, 0.5146, 0.4935, 0.4840, 0.4903,
% 0.5090, 0.5380, 0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891 ;
%
% Kdw() =  0.0510, 0.0405, 0.0331, 0.0278, 0.0242, 0.0217,
% 0.0200, 0.0189, 0.0182, 0.0178, 0.0176, 0.0176, 0.0179, 0.0193, 0.0224, 0.0280,
% 0.0369, 0.0498, 0.0526, 0.0577, 0.0640, 0.0723, 0.0842, 0.1065, 0.1578, 0.2409,
% 0.2892, 0.3124, 0.3296, 0.3290, 0.3559, 0.4105, 0.4278, 0.4521, 0.5116, 0.6514 ;
%
% l0;
% l1;
% k0;
% k1;
% m0;
% m1;
% kdiff;
% mdiff;
% num;
% den;
% frac;
% Kdw_l;
% M_l;
% Kd;
%
% int ref;
% int i;
%
% %   // -- INTERPOLATE TO WAVELENGTH OF INTEREST --  //
%
% for (i = 1; i < 36; i++)
%     if ( wave(i) >= lambda )
%         l1 = wave(i);
%         k1 = Kdw(i);
%         m1 = M(i);
%         l0 = wave(i-1);
%         k0 = Kdw(i-1);
%         m0 = M(i-1);
%         break;
%
%
%
%         num = lambda - l0;
%         den = l1 - l0;
%         frac = num / den;
%
%         kdiff = k1 - k0;
%         Kdw_l = k0 + frac*kdiff;
%
%         mdiff = m1 - m0;
%         M_l = m0 + frac*mdiff;
%
%
%         %   // -- GET REFERENCE WAVELENGTH (=490 FOR NOW) AND APPLY MODEL -- //
%
%         ref = 14;
%
%         Kd = (M_l/M(ref)) * (K490 - Kdw(ref)) + Kdw_l;
%
%         return Kd;
%



