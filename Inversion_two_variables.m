%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Combined_Claibration_BleachingModel.m
% Version     : 26.03.2019
% Input       : Known exposure age of the samples [year]
%               Experimental luminescence signal (Lx/Tx normalized)
%               Depth related to Exp. Lumi. Signal [mm]
% Output      : Bleaching rate SP [s-1]
%               Attenuation coeff. mu [mm-1]
% Inversion   : L1-norm waited over the experimental noise of the luminescence plateau
%               Inner loop is parallelized
% Coder       : B.Lehmann (lehmann.benj@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
addpath(genpath('/Users/blehman1/Documents/0_Pro/3_Publications/6_Application_OSL_Be_2_erosion/1_Traitement'));

TT   = 100;                   % Size of the inversion matrix(TT,TT)

t8   = 11*365.25*24.0*3600.0;   % Input the known age of the sample in years and conversion in seconds

%% Loading data file excel MBMV_.xlsx 

[num8,txt8,tab8]    =  xlsread('OSL_MBTP8_cal.xlsx');
n8                  =  length(num8(:,1));

%% Ordering of the imported data

n8C1              = 23;               % Disc number for core 1
n8C2              = 48;               % Disc number for core 2
x8(1:n8)          = (num8(1:n8,1));
[x_s8,ind8]       = sort(x8(:));

IR50_8(1:n8)      = (num8(1:n8,2));
IR50_e8(1:n8)     = (num8(1:n8,3));

IR50_s8           = IR50_8(ind8);
IR50_a8           = std(IR50_s8(47:n8));

%% Definition of parameter and domain of freedom 

x_model      = 0:0.25:25;

log_SP_max   = -3; 
log_SP_min   = -10;     
SP_max       = 10^(log_SP_max);   % s-1
SP_min       = 10^(log_SP_min);   % s-1
         
mu_max       = 4.0;               % mm-1
mu_min       = 0.1;               % mm-1

%% Initialization of the variable and matrix for the inversion

M8_IR50         = nan(TT,TT);

mu_matrix_IR50  = nan(TT,TT);
SP_matrix_IR50  = nan(TT,TT);

rand_vec1       = rand(TT,1);
rand_vec2       = rand(TT,1);

r_SP_IR50       = sort(10.^(log_SP_min+(log_SP_max-log_SP_min)*rand_vec1));
r_mu_IR50       = sort(mu_min+(mu_max-mu_min)*rand_vec2); 

%% Inversion

h               = waitbar(0,'Less than one song tato...');

    for i = 1:TT
        parfor j = 1:TT
            
            M8_IR50(i,j)   = 0;
                                
            IR50_8_th     = exp(-r_SP_IR50(j)*t8*exp(-r_mu_IR50(i)*x8));
            M8_IR50(i,j)  = sum((abs(IR50_8-IR50_8_th))./IR50_a8);                                   
                                    
            mu_matrix_IR50(i,j)  = r_mu_IR50(i);
            SP_matrix_IR50(i,j)  = r_SP_IR50(j);
            
        end
        waitbar(i/TT,h)
    end

close(h)

%% Transformation of the misfit M into likelihood chi and normalization

chi8_IR50       = 1./exp(0.5*M8_IR50);

max_chi8_IR50   = max(chi8_IR50(:));
norm_chi8_IR50  = chi8_IR50./max_chi8_IR50;

thr      = 0.01;

s_chi8_IR50        = norm_chi8_IR50(norm_chi8_IR50>thr);
s_SP8_IR50         = SP_matrix_IR50(norm_chi8_IR50>thr);
s_mu8_IR50         = mu_matrix_IR50(norm_chi8_IR50>thr);

s_chi8_BF        = norm_chi8_IR50(norm_chi8_IR50>0.95);
s_SP8_BF         = SP_matrix_IR50(norm_chi8_IR50>0.95);
s_mu8_BF         = mu_matrix_IR50(norm_chi8_IR50>0.95);
%% Extract 1d PDFs and confidence intervals for MBTP8

nbin               = 20;
% For SP
[n_SP_8,xout_SP_8]       =   hist(log10(s_SP8_IR50),nbin);
xwork_SP_8               =   cumsum(n_SP_8/sum(n_SP_8));

% -1 sigma 
ix_SP_8                  =   find(xwork_SP_8>0.175,1);
SP_1sd_8                 =   xout_SP_8(ix_SP_8);

% +1 sigma 
ix_SP_8                  =   find(xwork_SP_8>0.825,1);
SP_1su_8                 =   xout_SP_8(ix_SP_8);

% -2 sigma 
ix_SP_8                  =   find(xwork_SP_8>0.025,1);
SP_2sd_8                 =   xout_SP_8(ix_SP_8);

% +2 sigma 
ix_SP_8                  =   find(xwork_SP_8>0.925,1);
SP_2su_8                 =   xout_SP_8(ix_SP_8);

% Median
ix_SP_8                  =   find(xwork_SP_8>0.50,1);
SP_M_8                   =   xout_SP_8(ix_SP_8);

% Best fit
[n_SP_8_BF,xout_SP_8_BF] =   hist(log10(s_SP8_BF),nbin);
[xmin_SP_8_BF,tmin_SP_8] =   find(n_SP_8_BF == max(n_SP_8_BF),1);
SP_BF_8                  =   xout_SP_8_BF(tmin_SP_8);

%% For mu

[n_mu_8,xout_mu_8]       =   hist(s_mu8_IR50,nbin);
xwork_mu_8               =   cumsum(n_mu_8/sum(n_mu_8));

% -1 sigma 
ix_mu_8                  =   find(xwork_mu_8>0.175,1);
mu_1sd_8                 =   xout_mu_8(ix_mu_8);

% +1 sigma 
ix_mu_8                  =   find(xwork_mu_8>0.825,1);
mu_1su_8                 =   xout_mu_8(ix_mu_8);

% +2 sigma 
ix_mu_8                  =   find(xwork_mu_8>0.025,1);
mu_2sd_8                 =   xout_mu_8(ix_mu_8);

% +2 sigma 
ix_mu_8                  =   find(xwork_mu_8>0.925,1);
mu_2su_8                 =   xout_mu_8(ix_mu_8);

% Median
ix_mu_8                  =   find(xwork_mu_8>0.50,1);
mu_M_8                   =   xout_mu_8(ix_mu_8);

% Best fit
[n_mu_8_BF,xout_mu_8_BF]     =   hist(s_mu8_BF,nbin);
[xmin_mu_8_BF,tmin_mu_8_BF]  =   find(n_mu_8_BF == max(n_mu_8_BF),1);
mu_BF_8                      =   xout_mu_8_BF(tmin_mu_8_BF);

%% Print the results

fprintf('\nResult for the calibration with only MBTP8 \n')

SP_M_8_ok   = 10^(SP_M_8);
SP_BF_8_ok  = 10^(SP_BF_8);
SP_1su_8_ok = 10^(SP_1su_8);
SP_1sd_8_ok = 10^(SP_1sd_8);
SP_2su_8_ok = 10^(SP_2su_8);
SP_2sd_8_ok = 10^(SP_2sd_8);

disp(['SP Median     = ' num2str(SP_M_8_ok,3) ' s-1']);
disp(['SP BestFit    = ' num2str(SP_BF_8_ok,3) ' s-1']);
disp(['SP 1sigma sup = ' num2str(SP_1su_8_ok,3) ' s-1']);
disp(['SP 1sigma inf = ' num2str(SP_1sd_8_ok,3) ' s-1']);
disp(['SP 2sigma sup = ' num2str(SP_2su_8_ok,3) ' s-1']);
disp(['SP 2sigma inf = ' num2str(SP_2sd_8_ok,3) ' s-1']);

disp(['mu Median     = ' num2str(mu_M_8,3) ' mm-1']);
disp(['mu BestFit    = ' num2str(mu_BF_8,3) ' mm-1']);
disp(['mu 1sigma sup = ' num2str(mu_1su_8,3) ' mm-1']);
disp(['mu 1sigma inf = ' num2str(mu_1sd_8,3) ' mm-1']);
disp(['mu 2sigma sup = ' num2str(mu_2su_8,3) ' mm-1']);
disp(['mu 2sigma inf = ' num2str(mu_2sd_8,3) ' mm-1']);


%% Creating output models

IR50_8_M   = exp(-SP_M_8_ok*t8*exp(-mu_M_8*x_model));
IR50_8_BF  = exp(-SP_BF_8_ok*t8*exp(-mu_BF_8*x_model));

IR50_8_1sd_mu = exp(-SP_M_8_ok*t8*exp(-mu_1sd_8*x_model));
IR50_8_1su_mu = exp(-SP_M_8_ok*t8*exp(-mu_1su_8*x_model));
IR50_8_2sd_mu = exp(-SP_M_8_ok*t8*exp(-mu_2sd_8*x_model));
IR50_8_2su_mu = exp(-SP_M_8_ok*t8*exp(-mu_2su_8*x_model));

IR50_8_1sd_SP = exp(-SP_1sd_8_ok*t8*exp(-mu_M_8*x_model));
IR50_8_1su_SP = exp(-SP_1su_8_ok*t8*exp(-mu_M_8*x_model));
IR50_8_2sd_SP = exp(-SP_2sd_8_ok*t8*exp(-mu_M_8*x_model));
IR50_8_2su_SP = exp(-SP_2su_8_ok*t8*exp(-mu_M_8*x_model));


%% Plotting

figure('NumberTitle','off','Position',[00 00 1200 350])
sl =5;

subplot(1,2,1)
errorbar(x8(1:n8C1)     ,IR50_8(1:n8C1)     ,IR50_e8(1:n8C1)     ,'go' ,'markersize', sl)
hold on
errorbar(x8(n8C1+1:n8C2),IR50_8(n8C1+1:n8C2),IR50_e8(n8C1+1:n8C2),'b>' ,'markersize', sl)
errorbar(x8(n8C2+1:n8)  ,IR50_8(n8C2+1:n8)  ,IR50_e8(n8C2+1:n8)  ,'rd' ,'markersize', sl)
plot(x_model,IR50_8_M ,'k', 'LineWidth',1)
plot(x_model,IR50_8_BF,'k--', 'LineWidth',2)
plot(x_model,IR50_8_1su_mu ,'Color', [0.7 0.7 0.7], 'LineWidth',2)
plot(x_model,IR50_8_1sd_mu ,'Color', [0.7 0.7 0.7], 'LineWidth',2)
plot(x_model,IR50_8_1su_SP ,'Color', [0.9 0.9 0.9], 'LineWidth',2)
plot(x_model,IR50_8_1sd_SP ,'Color', [0.9 0.9 0.9], 'LineWidth',2)

%plot(x_model,IR50_8_2su_mu ,'Color', [0.7 0.7 0.7], 'LineWidth',2)
%plot(x_model,IR50_8_2sd_mu ,'Color', [0.7 0.7 0.7], 'LineWidth',2)


xlabel('Depth [mm]')
ylabel('IRSL Normalized Intensity')
legend('Core 1','Core 2','Core 3','MBTP8 Median','MBTP8 BestFit','± 1sigma of mu','± 1sigma of SP','Location','Southeast')
title('MBTP8 - IR50')
axis([0 28 0 1.3])

subplot(1,2,2)
surface(r_SP_IR50,r_mu_IR50,norm_chi8_IR50); axis square; colorbar;shading interp;set(gca,'XScale','log');
hold on

SP_vec_plot      = 0:0.01:4;
SP_M_8_line_plot = SP_M_8_ok.*ones(size(SP_vec_plot));
SP_BF_8_line_plot = SP_BF_8_ok.*ones(size(SP_vec_plot));

mu_vec_plot      = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3];
mu_M_8_line_plot = mu_M_8.*ones(size(mu_vec_plot));
mu_BF_8_line_plot = mu_BF_8.*ones(size(mu_vec_plot));

plot(mu_vec_plot,mu_M_8_line_plot,'w','LineWidth',1)
plot(SP_M_8_line_plot,SP_vec_plot,'w','LineWidth',1)

plot(mu_vec_plot,mu_BF_8_line_plot,'w--','LineWidth',1)
plot(SP_BF_8_line_plot,SP_vec_plot,'w--','LineWidth',1)

xlabel('Bleaching rate [s-1]')
ylabel('Attenuation coeff. [mm-1]')
axis([1e-10 1e-3 0.1 4])
title('Cal. MBTP8- IR50')

