%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : B_Inversion_Erosion_ts_MBTP1.m
%
% Version     : 17.07.2019
%
% Coder       : Benjamin Lehmann (lehmann.benj@gmail.com)
%
% Aim         : Invert pair of ts and e_r which could solve TCN and OSL system
%               From 10Be concentration and OSL signal
%               and considering erosion histories evolving as step function (ts: erosion time; e_r: erosion rate)
%
% Input       : x_data     = Depth at which luminescence signal is measure [mm] (is coverted into m later)
%               L_data     = Luminescence signal
%               L_er_data  = Error on luminescence signal
%               Ddot_input = Dose rate [Gy/a]
%               sp0        = Bleaching rate [s-1]
%               mu0        = Attenuation coeff. mu [mm-1]
%               D0         = Environmental dose [Gy] 
%               All variable and matrix from A_Mask_cosmo.m
%
% Output      : norm_chi    = Probability distribution for ts and e_r 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc;  
close all;

% Set matlab passes
% addpath('...');

tic()

SampleName     = 'MBTP1';

%% Inputs for OSL
% Experiemental values
OSL_SampleName = ['OSL_IR50_' SampleName];
[num]          = xlsread(OSL_SampleName);
x_data         = num(:,1)*1e-3;           % data in [mm] transform in [m]
L_data         = num(:,2);
L_er_data      = num(:,3);
Ddot_input     = num(1,5);

% Calibration values
sp0            = 4.1E-06;                 % [s-1]  from calibration on samples MBTP8 & MBTP6
mu0            = 0.596;                   % [mm-1] from calibration on samples MBTP8 & MBTP6
sp             = sp0*3600*24*365;         % [a-1]
mu             = mu0*1e3;                 % [m-1]
Ddot           = Ddot_input/1e3;          % [Gy/a]
D0             = 500;                     % [Gy] from literature, could be measured for each sample

%% Load mask

load([SampleName '_TCN_Mask.mat']);

% Mask_cosmo_corr  = Mask_cosmo;
% t_corr_mat_corr  = t_corr_mat;
% 
% for j=1:55 
%     for i=70:100
%         if Mask_cosmo_corr(i,j)==0
%         Mask_cosmo_corr(i,j)   = 1;
%         end        
%     end
% end
% 
% for j=1:55
%     for i=70:100
%      if isnan(t_corr_mat_corr(i,j))==1
%      t_corr_mat_corr(i,j) = min(min(t_corr_mat));
%      end
%     end
% end


%% Numerics: Discretization, time steps and plotting frequency

dt              = 0.01;                           %  Time step [a]
disp(['dt = ' num2str(dt) ' a']);

[xs_data,ind_x] = sort(x_data(:));
Ls_data         = L_data(ind_x);
ind_a           = min(find(Ls_data>0.9));         % Find the index when the plateau starts 
a               = std(Ls_data(ind_a:end));        % Calculate the std over the plateau

xmax            = max(x_data);                    % [m]
n               = 101;
dx              = xmax/(n-1);
x_final         = 0:dx:xmax;
L_int           = interp1(x_data,L_data,x_final);

M               = nan(TT,TT);                     % Misfit matrix

%% Inversion calculation

disp('The inversion strats NOW')
% For Cluster work
% pool = parpool(20);

for i = 1:TT                                            % loop over e_rpar
disp(['Outer loop # ' num2str(i)])
parfor j = 1:TT                                         % loop over tau_r
    
%    if Mask_cosmo_corr(i,j) == 0
    if Mask_cosmo(i,j) == 0
       M(i,j)          = NaN;

    else
%    t                = t_corr_mat_corr(i,j);
    t                = t_corr_mat(i,j);
    nt               = floor(t/dt);
    xth              = x_final;
    Lth              = ones(1,n);
    Lth(1,1)         = 0;
    time             = zeros(nt,1);
    time(1)          = 0;
    M(i,j)           = 0;
    Lthp             = Lth;
    res              = 0;
    
    for it           = 2:nt
        time(it)     = time(it-1)+dt;

        if time(it) <= t-ts_r(j) 
        e_rate       = 0; 
        else
        e_rate       = e_r(i);
        end
        
        advection      = e_rate*(-Lth(3:end)+4*Lth(2:end-1)-3*Lth(1:end-2))/(2*dx);
        Lth(1:n-2)     = Lth(1:n-2)+dt.*(Ddot/D0*(1-Lth(1:n-2))-(Lth(1:n-2)*sp).*exp(-mu.*xth(1:n-2)))+advection*dt;  % Eq. 1 Lehmann et al., 2019
        Lth            = max(Lth,0);
        Lth(end-1:end) = 1; 
        Lth(end)       = 1; 
            
        res = sum(abs(Lthp-Lth));
        if (res<1e-8)
        break
        end
        Lthp=Lth;
    end

    M(i,j) = nansum((Lth-L_int).^2/a^2);  % Square-integrable function L2 norm            

   end
end
end
disp('The Inversion is OVER') 

%% Transformation of the misfit into likelihod and normalization

chi       = 1./exp(0.5*M) ;
max_chi   = max(max(chi)) ;
norm_chi  = chi./max_chi  ;

%% Saving

name_mat = [SampleName '_inversion_result.mat'];
save(name_mat);
disp('Saving successfully proceed');
toc()

% If pool where use --> necessity to close it
% delete(pool)
