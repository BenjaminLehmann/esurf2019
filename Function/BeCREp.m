% Function BeCREp 
% Inputs:
%     SampleName
%     Erosion_ct_on_off : if 0 Erosion is constant 
%                         if 1 Erosion is following a step function in time
%                         (er_max during ts)
%     erosion : erosion rate [m/a] is Erosion_ct_on_off = 0
%     er_max  : erosion rate [m/a] is Erosion_ct_on_off = 1
%     ts : time [a] at which er_max start to be 
% Outputs
%     Age  : 10Be exposure age [a]
%     Err  : internal error on Age
%     Err2 : external error on Age
%     Erosion_rate : erosion rate using 10Be conc. assuming erosion steady
%     state Eq. 11 Lal et al., 1991)
%     Convergence : if 0 = no convergence, if 1 = convergence


function [Age,Err,Err2,Erosion_rate,Convergence]=BeCREp(SampleName,Erosion_ct_on_off,erosion,er_max,ts)

[num]       = xlsread(SampleName);
Convergence = NaN;
%% INPUT caracteristics of the sample

VecLat     = num(:,1);           % Laltitude of the sample (+) to the north (-) to the south [deg]
VecLon     = num(:,2);           % Longitude of the sample (+) to the est (-) to the west [deg]
VecAlt     = num(:,3);           % Altitude of the sample [m.a.s.l.]
VecConc    = num(:,4);           % Nuclide concentration of the sample [at/g]
VecErrConc = num(:,5);           % 1 sigma of the nuclide concentration of the sample [at/g]
VecShield  = num(:,6);           % Shielding factor of the sample
VecDens    = num(:,7)*1e6;       % Density of the sample [g/cm3] + conversion into [g/m3]
VecThick   = num(:,8)/100;       % Thickness of the sample [cm] + conversion into [m]
VecErosion = erosion/100 ;       % Erosion [cm/a] + conversion into [m/a]

%% DO NOT CHANGE BELOW
% Corrections the density, thickness and shielding

Attlg      = 160*1e4;                            % Attenuation length for spallation in rock [g/cm2] (Gosse and Phillips, 2001) + conversion into [g/m2]
ThickCorr  = (Attlg/(VecThick*VecDens))*(1-exp(-1*(VecThick*VecDens)/Attlg));  
CorrConc   = VecConc/(ThickCorr*VecShield);      % without snow correction
VecMu      = VecDens/Attlg;

% Caracteristics of the local production rate (PR)

LocalLatV  = 46.42083;         % Laltitude of the Chironico landslide [deg]
LocalLonV  = 8.850456;         % Longitude of the Chironico landslide [deg]
LocalAltV  = 786.6667;         % Altitude of the Chironico landslide [m.a.s.l.]
LocalAgeV  = 13378;            % Age of the Chironico landslide [a]
SelPR(1)   = 4.1572;           % Chironico landslide local PR [at/g/yr]
SelPR(2)   = 0.1003;           % 1sigma on the local PR [at/g/yr]
tBe        = 1387000;          % 10Be half-life [a] (Nishiizumi et al., 2007)
Lambda10Be = log(2)/tBe;

% Flags for BLARD routines
Atm        = 0;                % O stand for the ERA40 Atmospheric model
Nucl       = 10;               % 10 stand for the 10Be nuclide 

% Age calculation using BLARD routines
load('GMDB.mat')
NumGMDB      = GMDB.GLOPIS;
[VecT,VecSF] = LSDv9(LocalLatV,LocalLonV,LocalAltV,Atm,LocalAgeV,-1,10,NumGMDB);
[~,SF]       = LSDv9(VecLat,VecLon,VecAlt,Atm,0,-1,Nucl,NumGMDB);

%% Choose if erosion is constant in time or not 
% Erosion_ct_on_off== 0 : constant
% Erosion_ct_on_off== 1 : non-constant

P0         = SelPR(1)*SF;
P0_1sigma  = SelPR(2)*SF;

if Erosion_ct_on_off== 0   
    Convergence = 1;
    CosmoAge    = (-1/(Lambda10Be+VecMu*VecErosion))*log(1-((Lambda10Be+VecMu*VecErosion)*CorrConc/(SelPR(1)*SF)))/1000;
    else

    % Find new age with erosion function
   
    CosmoAge   = 20000;         % Initial guess [a]
    residuals  = 1e6;
    residualsp = residuals+1;
    nt         = 100;
    lambda     = (1.5102/P0);   % lambda should <= 0.1
%     disp(lambda)
    for it    = 1:nt
        tendp           = CosmoAge;
        %fprintf('%d %d %d %d \n',CosmoAge,P0,er_max,ts)
        [concentration] = cosmo_conc2age(CosmoAge,P0,er_max,ts,Lambda10Be,VecMu);
        CosmoAge        = CosmoAge-lambda*(concentration-CorrConc);
        residuals       = abs(CosmoAge-tendp);

        figure(4)
        plot(it,abs(CosmoAge-tendp),'o');
        hold on
        xlabel('Iteration')
        ylabel('Residual [Years]')
        drawnow

        if (residuals < 1.)
            Convergence   = 1;
          %   disp(it)
            break
        end
        
        if (residuals > 1.) && (it == nt) 
            %disp('this is not going to work!')
            Convergence   = 0;
            break
        end
        
        % Optimisation settings
        
        if (residuals > 1e03) && (it == nt/5) 
            %disp('this is not going to work!')
            Convergence   = 0;
            break
        end
        
        if (residuals <= residualsp)                  
            residualsp = residuals;
        elseif residuals-residualsp > 100
            Convergence   = 0;
            break
        end
    end
    CosmoAge   = CosmoAge/1000;

    if  isreal(CosmoAge) == 0
        ErrPRflag        =  1;
        CosmoAge         = -100;
    end
end

%disp(P0)
%disp(P0_1sigma)

if Convergence     == 0
    Age                = NaN;      
    Err                = NaN;
    Err2               = NaN;
    Erosion_rate       = (SelPR(1)*SF*exp(-VecMu*VecThick)-Lambda10Be*CorrConc)/(CorrConc*VecMu);    % Eq. 11 Lal et al., 1991
else
    ErrCosmo           = CosmoAge*VecErrConc/VecConc;
    [Age,Err,Err2,X,Y] = AgeCosmoAgeReelLSD(CosmoAge,ErrCosmo,(SelPR(2)/SelPR(1)),VecLat,VecLon,VecAlt,NumGMDB,Atm,Nucl);
    Age                = Age*1e3;       %  conversion ka to a
    Err                = Err*1e3;
    Err2               = Err2*1e3;
    Erosion_rate       = (SelPR(1)*SF*exp(-VecMu*VecThick)-Lambda10Be*CorrConc)/(CorrConc*VecMu);   % Eq. 11 Lal et al., 1991
    Convergence        = 1;   
end

end
