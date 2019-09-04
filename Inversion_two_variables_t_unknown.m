%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Inversion_two_variables_t_unknown.m
% Version     : 27.08.2019
% Input       : Unknown exposure age of the samples [year]
%               Experimental luminescence signal (Lx/Tx normalized)
%               Depth related to Exp. Lumi. Signal [mm]
% Output      : Bleaching rate SP x exposure time 
%               Attenuation coeff. mu [mm-1]
% Inversion   : L1-norm waited over the experimental noise of the luminescence plateau
% Contact     : B.Lehmann (lehmann.benj@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
addpath(genpath('/Users/blehman1/Documents/0_Pro/5_Code'));  % Set the depository path to the actual working folder

TT      = 1000;                                              % Size of the inversion matrix(TT,TT)

%% Loading data file excel MBMV_.xlsx 

[num,txt,tab]    =  xlsread('OSL_MBTP7_cal.xlsx');
n                =  length(num(:,1));

%% Ordering of the imported data

nC1          = 23;                     % Disc number for core 1
nC2          = 47;                     % Disc number for core 2
ind_plt      = 64;                     % at which index the plateaus starts

x(1:n)       = (num(1:n,1));           % 1st colomn of excell table = depth [mm]
IR50(1:n)    = (num(1:n,2));           % 2nd colomn of excell table = Lx/Tx
IR50_e(1:n)  = (num(1:n,3));           % 3rd colomn of excell table = Error Lx/Tx

[x_s,ind]    = sort(x(:));             % sorting increasingly the depth 
IR50_s       = IR50(ind);              % sorting Lx/Tx with increasing depth
IR50_a       = std(IR50_s(ind_plt:n)); % calculating the standard deviation on the plateau

%% Definition of parameter and domain of freedom 

x_model      = 0:0.25:25;

log_SP_t_max   = 4; 
log_SP_t_min   = 1;     
SP_t_max       = 10^(log_SP_t_max);   % Dimensionless
SP_t_min       = 10^(log_SP_t_min);   % Dimensionless
         
mu_max       = 2.0;               % mm-1
mu_min       = 0.1;               % mm-1

%% Initialization of the variable and matrix for the inversion

M            = nan(TT,TT);    % creation of an empty misfit matrix

mu_matrix    = nan(TT,TT);    % creation of an empty mu matrix
SP_t_matrix  = nan(TT,TT);    % creation of an empty SP_t matrix

rand_vec1    = rand(TT,1);    % random vector of size TT to sample SP_t
rand_vec2    = rand(TT,1);    % random vector of size TT to sample mu

r_SP_t       = sort(10.^(log_SP_t_min+(log_SP_t_max-log_SP_t_min)*rand_vec1));    % r_SP_t is randomly sampled in log space
r_mu         = sort(mu_min+(mu_max-mu_min)*rand_vec2);                            % mu is randomly sampled in normal space

%% Inversion

h               = waitbar(0,'Put a good song and relax');  % creation of the loading bar (just for fun)

    for i = 1:TT
        for j = 1:TT                                       % "parfor" for paralellelization, for normal claculation put "for"
            
            M(i,j)   = 0;                                  % initialization of the misfit at i and j index
                                
            IR50_th  = exp(-r_SP_t(j)*exp(-r_mu(i)*x));    % creation of synthetic luminescence signal
            M(i,j)   = sum((abs(IR50-IR50_th))./IR50_a);   % calculation of the misfit between experimental data and synthetic signal                             
                                    
            mu_matrix(i,j)  = r_mu(i);                     % compile mu tested into matrix
            SP_t_matrix(i,j)  = r_SP_t(j);                 % compile SP_t tested into matrix
            
        end
        waitbar(i/TT,h)                                    % incrementation of the loading bar (just for fun)
    end

close(h)                                                   % close the loading bar

%% Transformation of the misfit M into likelihood chi and normalization

chi        = 1./exp(0.5*M);

max_chi    = max(chi(:));
norm_chi   = chi./max_chi;

thr        = 0.01;

s_chi      = norm_chi(norm_chi>thr);
s_SP_t     = SP_t_matrix(norm_chi>thr);
s_mu       = mu_matrix(norm_chi>thr);

s_chi_BF   = norm_chi(norm_chi>0.95);
s_SP_t_BF  = SP_t_matrix(norm_chi>0.95);
s_mu_BF    = mu_matrix(norm_chi>0.95);

%% Extract 1d PDFs and confidence intervals for MBTP8

nbin                 = 20;
% For SP
[n_SP_t,xout_SP_t]   =   hist(log10(s_SP_t),nbin);
xwork_SP_t           =   cumsum(n_SP_t/sum(n_SP_t));

% -1 sigma 
ix_SP_t              =   find(xwork_SP_t>0.175,1);
SP_t_1sd               =   xout_SP_t(ix_SP_t);

% +1 sigma 
ix_SP_t                =   find(xwork_SP_t>0.825,1);
SP_t_1su               =   xout_SP_t(ix_SP_t);

% -2 sigma 
ix_SP_t                =   find(xwork_SP_t>0.025,1);
SP_t_2sd               =   xout_SP_t(ix_SP_t);

% +2 sigma 
ix_SP_t                =   find(xwork_SP_t>0.925,1);
SP_t_2su               =   xout_SP_t(ix_SP_t);

% Median
ix_SP_t                =   find(xwork_SP_t>0.50,1);
SP_t_M                 =   xout_SP_t(ix_SP_t);

% Best fit
[n_SP_t_BF,xout_SP_BF] =   hist(log10(s_SP_t_BF),nbin);
[xmin_SP_t_BF,tmin_SP_t] =   find(n_SP_t_BF == max(n_SP_t_BF),1);
SP_t_BF                =   xout_SP_BF(tmin_SP_t);

%% For mu

[n_mu,xout_mu]           =   hist(s_mu,nbin);
xwork_mu                 =   cumsum(n_mu/sum(n_mu));

% -1 sigma 
ix_mu                    =   find(xwork_mu>0.175,1);
mu_1sd                   =   xout_mu(ix_mu);

% +1 sigma 
ix_mu                    =   find(xwork_mu>0.825,1);
mu_1su                   =   xout_mu(ix_mu);

% +2 sigma 
ix_mu                    =   find(xwork_mu>0.025,1);
mu_2sd                   =   xout_mu(ix_mu);

% +2 sigma 
ix_mu                    =   find(xwork_mu>0.925,1);
mu_2su                   =   xout_mu(ix_mu);

% Median
ix_mu                    =   find(xwork_mu>0.50,1);
mu_M                     =   xout_mu(ix_mu);

% Best fit
[n_mu_BF,xout_mu_BF]     =   hist(s_mu_BF,nbin);
[xmin_mu_BF,tmin_mu_BF]  =   find(n_mu_BF == max(n_mu_BF),1);
mu_BF                    =   xout_mu_BF(tmin_mu_BF);

%% Print the results

fprintf('\nResult for the calibration with only MBTP8 \n')

SP_t_M_ok   = 10^(SP_t_M);
SP_t_BF_ok  = 10^(SP_t_BF);
SP_t_1su_ok = 10^(SP_t_1su);
SP_t_1sd_ok = 10^(SP_t_1sd);
SP_t_2su_ok = 10^(SP_t_2su);
SP_t_2sd_ok = 10^(SP_t_2sd);

disp(['SPxt Median     = ' num2str(SP_t_M_ok,3) ]);
disp(['SPxt BestFit    = ' num2str(SP_t_BF_ok,3)]);
disp(['SPxt 1sigma sup = ' num2str(SP_t_1su_ok,3)]);
disp(['SPxt 1sigma inf = ' num2str(SP_t_1sd_ok,3)]);
disp(['SPxt 2sigma sup = ' num2str(SP_t_2su_ok,3)]);
disp(['SPxt 2sigma inf = ' num2str(SP_t_2sd_ok,3)]);

disp(['mu Median     = ' num2str(mu_M,3) ' mm-1']);
disp(['mu BestFit    = ' num2str(mu_BF,3) ' mm-1']);
disp(['mu 1sigma sup = ' num2str(mu_1su,3) ' mm-1']);
disp(['mu 1sigma inf = ' num2str(mu_1sd,3) ' mm-1']);
disp(['mu 2sigma sup = ' num2str(mu_2su,3) ' mm-1']);
disp(['mu 2sigma inf = ' num2str(mu_2sd,3) ' mm-1']);


%% Creating output models

IR50_M   = exp(-SP_t_M_ok*exp(-mu_M*x_model));
IR50_BF  = exp(-SP_t_BF_ok*exp(-mu_BF*x_model));

%% Plotting

figure('NumberTitle','off','Position',[00 00 1200 350])
sl =5;

subplot(1,2,1)
errorbar(x(1:nC1)     ,IR50(1:nC1)     ,IR50_e(1:nC1)     ,'go' ,'markersize', sl)
hold on
errorbar(x(nC1+1:nC2),IR50(nC1+1:nC2),IR50_e(nC1+1:nC2),'b>' ,'markersize', sl)
errorbar(x(nC2+1:n)  ,IR50(nC2+1:n)  ,IR50_e(nC2+1:n)  ,'rd' ,'markersize', sl)
plot(x_model,IR50_M ,'k', 'LineWidth',1)
plot(x_model,IR50_BF,'k--', 'LineWidth',2)

xlabel('Depth [mm]')
ylabel('IRSL Normalized Intensity')
legend('Core 1','Core 2','Core 3','MBTP7 Median','MBTP7 BestFit','Location','Southeast')
title('MBTP7 - IR50')
axis([0 28 0 1.3])

subplot(1,2,2)
surface(r_SP_t,r_mu,norm_chi); axis square; colorbar;shading interp;set(gca,'XScale','log');
hold on

SP_vec_plot      = 0:0.01:4;
SP_M_line_plot   = SP_t_M_ok.*ones(size(SP_vec_plot));
SP_BF_line_plot  = SP_t_BF_ok.*ones(size(SP_vec_plot));

mu_vec_plot      = [1e1,1e2,1e3,1e4];
mu_M_line_plot   = mu_M.*ones(size(mu_vec_plot));
mu_BF_line_plot  = mu_BF.*ones(size(mu_vec_plot));

plot(mu_vec_plot,mu_M_line_plot,'w','LineWidth',1)
plot(SP_M_line_plot,SP_vec_plot,'w','LineWidth',1)

plot(mu_vec_plot,mu_BF_line_plot,'w--','LineWidth',1)
plot(SP_BF_line_plot,SP_vec_plot,'w--','LineWidth',1)

xlabel('Bleaching rate x time [Dimensionless]')
ylabel('Attenuation coeff. [mm-1]')
axis([SP_t_min SP_t_max mu_min mu_max])
title('Cal. MBTP7- IR50')

