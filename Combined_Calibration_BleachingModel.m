%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Combined_Claibration_BleachingModel.m
% Version     : 27.08.2019
% Input       : Known exposure age of the samples [year]
%               Experimental luminescence signal (Lx/Tx normalized)
%               Depth related to Exp. Lumi. Signal [mm]
% Output      : Bleaching rate SP [s-1]
%               Attenuation coeff. mu [mm-1]
% Inversion   : L1-norm waited over the experimental noise of the luminescence plateau
%               Inner loop is parallelized
% Contact     : B.Lehmann (lehmann.benj@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
addpath(genpath('/Users/blehman1/Documents/0_Pro/3_Publications/6_Application_OSL_Be_2_erosion/1_Traitement'));

TT   = 1000;                    % Size of the inversion matrix(TT,TT)

t7   = 2*365.25*24.0*3600.0;    % Input the known age of the sample in years and conversion in seconds
t8   = 11*365.25*24.0*3600.0;   % Input the known age of the sample in years and conversion in seconds

%% Loading data file excel MBMV_.xlsx 

[num7,txt7,tab7]    =  xlsread('OSL_MBTP7_cal.xlsx');
n7                  =  length(num7(:,1));

[num8,txt8,tab8]    =  xlsread('OSL_MBTP8_cal.xlsx');
n8                  =  length(num8(:,1));

%% Ordering of the imported data
% For sample MBTP7
n7C1              = 23;                         % Disc number for core 1
n7C2              = 47;                         % Disc number for core 2
x7(1:n7)          = (num7(1:n7,1));             % 1st colomn of excell table of sample MBTP7 = depth [mm]
IR50_7(1:n7)      = (num7(1:n7,2));             % 2nd colomn of excell table of sample MBTP7 = Lx/Tx 
IR50_e7(1:n7)     = (num7(1:n7,3));             % 3rd colomn of excell table of sample MBTP7 = Error Lx/Tx

[x_s7,ind7]       = sort(x7(:));                % sorting increasingly the depth 
IR50_s7           = IR50_7(ind7);               % sorting Lx/Tx with increasing depth
ind_plt_7         = 64;                         % at which index the plateaus starts
IR50_a7           = std(IR50_s7(ind_plt_7:n7)); % at which index the plateaus starts

% Same for sample MBTP8
n8C1              = 23;                         
n8C2              = 48;                         
x8(1:n8)          = (num8(1:n8,1));
IR50_8(1:n8)      = (num8(1:n8,2));
IR50_e8(1:n8)     = (num8(1:n8,3));

[x_s8,ind8]       = sort(x8(:));
IR50_s8           = IR50_8(ind8);
ind_plt_8         = 47;
IR50_a8           = std(IR50_s8(ind_plt_8:n8));

%% Definition of parameter and domain of freedom 

x_model      = 0:0.25:25;

log_SP_max   = -3; 
log_SP_min   = -10;     
SP_max       = 10^(log_SP_max);   % s-1
SP_min       = 10^(log_SP_min);   % s-1
         
mu_max       = 4.0;               % mm-1
mu_min       = 0.1;               % mm-1

%% Initialization of the variable and matrix for the inversion

M7         = nan(TT,TT);
M8         = nan(TT,TT);

mu_matrix  = nan(TT,TT);
SP_matrix  = nan(TT,TT);

rand_vec1  = rand(TT,1);
rand_vec2  = rand(TT,1);

r_SP       = sort(10.^(log_SP_min+(log_SP_max-log_SP_min)*rand_vec1));
r_mu       = sort(mu_min+(mu_max-mu_min)*rand_vec2); 

%% Inversion

h               = waitbar(0,'Less than one song tato...');

    for i = 1:TT
        parfor j = 1:TT
            
            M8(i,j)   = 0;
            M7(i,j)   = 0;
                                
            IR50_8_th     = exp(-r_SP(j)*t8*exp(-r_mu(i)*x8));
            M8(i,j)  = sum((abs(IR50_8-IR50_8_th))./IR50_a8);                                   
            
            IR50_7_th     = exp(-r_SP(j)*t7*exp(-r_mu(i)*x7));
            M7(i,j)  = sum((abs(IR50_7-IR50_7_th))./IR50_a7);                                   
                        
            mu_matrix(i,j)  = r_mu(i);
            SP_matrix(i,j)  = r_SP(j);
            
        end
        waitbar(i/TT,h)
    end

close(h)

%% Transformation of the misfit M into likelihood chi and normalization

chi7          = 1./exp(0.5*M7);    % Transformation misfit into likelihood
chi8          = 1./exp(0.5*M8);

max_chi7      = max(chi7(:));
norm_chi7     = chi7./max_chi7;

max_chi8      = max(chi8(:));
norm_chi8     = chi8./max_chi8;

chi_all       = norm_chi7.*norm_chi8;
max_chi_all   = max(chi_all(:));
norm_chi_all  = chi_all./max_chi_all;

thr           = 0.01;

s_chi7        = norm_chi7(norm_chi7>thr);
s_SP7         = SP_matrix(norm_chi7>thr);
s_mu7         = mu_matrix(norm_chi7>thr);

s_chi8        = norm_chi8(norm_chi8>thr);
s_SP8         = SP_matrix(norm_chi8>thr);
s_mu8         = mu_matrix(norm_chi8>thr);

s_chi_all     = norm_chi_all(norm_chi_all>thr);
s_SP_all      = SP_matrix(norm_chi_all>thr);
s_mu_all      = mu_matrix(norm_chi_all>thr);

s_chi7_BF     = norm_chi7(norm_chi7>0.95);
s_SP7_BF      = SP_matrix(norm_chi7>0.95);
s_mu7_BF      = mu_matrix(norm_chi7>0.95);

s_chi8_BF     = norm_chi8(norm_chi8>0.95);
s_SP8_BF      = SP_matrix(norm_chi8>0.95);
s_mu8_BF      = mu_matrix(norm_chi8>0.95);

s_chi_all_BF  = norm_chi_all(norm_chi_all>0.95);
s_SP_all_BF   = SP_matrix(norm_chi_all>0.95);
s_mu_all_BF   = mu_matrix(norm_chi_all>0.95);

%% Extract 1d PDFs and confidence intervals for combined signal
nbin               = 20;
% For SP
[n_SP_all,xout_SP_all]     =   hist(log10(s_SP_all),nbin);
xwork_SP_all               =   cumsum(n_SP_all/sum(n_SP_all));

ix_SP_all                  =   find(xwork_SP_all>0.175,1);
SP_1sd_all                 =   xout_SP_all(ix_SP_all);

ix_SP_all                  =   find(xwork_SP_all>0.825,1);
SP_1su_all                 =   xout_SP_all(ix_SP_all);

ix_SP_all                  =   find(xwork_SP_all>0.025,1);
SP_2sd_all                 =   xout_SP_all(ix_SP_all);

ix_SP_all                  =   find(xwork_SP_all>0.925,1);
SP_2su_all                 =   xout_SP_all(ix_SP_all);

ix_SP_all                  =   find(xwork_SP_all>0.50,1);
SP_M_all                   =   xout_SP_all(ix_SP_all);

% [xmin_SP_all,tmin_SP_all]  =   find(n_SP_all == max(n_SP_all),1);
% SP_BF_all                  =   xout_SP_all(tmin_SP_all);

[n_SP_BF_all,xout_SP_all_BF] =   hist(log10(s_SP_all_BF),nbin);
[xmin_SP_BF,tmin_SP_all_BF]  =   find(n_SP_BF_all == max(n_SP_BF_all),1);
SP_BF_all                    =   xout_SP_all_BF(tmin_SP_all_BF);


% For mu
[n_mu_all,xout_mu_all]     =   hist(s_mu_all,nbin);
xwork_mu_all               =   cumsum(n_mu_all/sum(n_mu_all));

ix_mu_all                  =   find(xwork_mu_all>0.175,1);
mu_1sd_all                 =   xout_mu_all(ix_mu_all);

ix_mu_all                  =   find(xwork_mu_all>0.825,1);
mu_1su_all                 =   xout_mu_all(ix_mu_all);

ix_mu_all                  =   find(xwork_mu_all>0.025,1);
mu_2sd_all                 =   xout_mu_all(ix_mu_all);

ix_mu_all                  =   find(xwork_mu_all>0.925,1);
mu_2su_all                 =   xout_mu_all(ix_mu_all);

ix_mu_all                  =   find(xwork_mu_all>0.50,1);
mu_M_all                   =   xout_mu_all(ix_mu_all);

[n_mu_BF_all,xout_mu_all_BF]   =   hist(s_mu_all_BF,nbin);
[xmin_mu_BF,tmin_mu_all_BF]    =   find(n_mu_BF_all == max(n_mu_BF_all),1);
mu_BF_all                      =   xout_mu_all_BF(tmin_mu_all_BF);

%% Extract 1d PDFs and confidence intervals for MBTP7
nbin               = 20;

% For SP
[n_SP_7,xout_SP_7]     =   hist(log10(s_SP7),nbin);
xwork_SP_7               =   cumsum(n_SP_7/sum(n_SP_7));

ix_SP_7                  =   find(xwork_SP_7>0.175,1);
SP_1sd_7                 =   xout_SP_7(ix_SP_7);

ix_SP_7                  =   find(xwork_SP_7>0.825,1);
SP_1su_7                 =   xout_SP_7(ix_SP_7);

ix_SP_7                  =   find(xwork_SP_7>0.025,1);
SP_2sd_7                 =   xout_SP_7(ix_SP_7);

ix_SP_7                  =   find(xwork_SP_7>0.925,1);
SP_2su_7                 =   xout_SP_7(ix_SP_7);

ix_SP_7                  =   find(xwork_SP_7>0.50,1);
SP_M_7                   =   xout_SP_7(ix_SP_7);

[n_SP_BF_7,xout_SP_7_BF]    =   hist(log10(s_SP7_BF),nbin);
[xmin_SP_BF_7,tmin_SP_7_BF] =   find(n_SP_BF_7 == max(n_SP_BF_7),1);
SP_BF_7                     =   xout_SP_7_BF(tmin_SP_7_BF);

% For mu
[n_mu_7,xout_mu_7]     =   hist(s_mu7,nbin);
xwork_mu_7               =   cumsum(n_mu_7/sum(n_mu_7));

ix_mu_7                  =   find(xwork_mu_7>0.175,1);
mu_1sd_7                 =   xout_mu_7(ix_mu_7);

ix_mu_7                  =   find(xwork_mu_7>0.825,1);
mu_1su_7                 =   xout_mu_7(ix_mu_7);

ix_mu_7                  =   find(xwork_mu_7>0.025,1);
mu_2sd_7                 =   xout_mu_7(ix_mu_7);

ix_mu_7                  =   find(xwork_mu_7>0.925,1);
mu_2su_7                 =   xout_mu_7(ix_mu_7);

ix_mu_7                  =   find(xwork_mu_7>0.50,1);
mu_M_7                   =   xout_mu_7(ix_mu_7);

[n_mu_BF_7,xout_mu_7_BF]  =   hist(s_mu7_BF,nbin);
[xmin_mu_BF,tmin_mu_7_BF] =   find(n_mu_BF_7 == max(n_mu_BF_7),1);
mu_BF_7                   =   xout_mu_7_BF(tmin_mu_7_BF);

%% Extract 1d PDFs and confidence intervals for MBTP8

nbin               = 20;
% For SP
[n_SP_8,xout_SP_8]       =   hist(log10(s_SP8),nbin);
xwork_SP_8               =   cumsum(n_SP_8/sum(n_SP_8));

ix_SP_8                  =   find(xwork_SP_8>0.175,1);
SP_1sd_8                 =   xout_SP_8(ix_SP_8);

ix_SP_8                  =   find(xwork_SP_8>0.825,1);
SP_1su_8                 =   xout_SP_8(ix_SP_8);

ix_SP_8                  =   find(xwork_SP_8>0.025,1);
SP_2sd_8                 =   xout_SP_8(ix_SP_8);

ix_SP_8                  =   find(xwork_SP_8>0.925,1);
SP_2su_8                 =   xout_SP_8(ix_SP_8);

ix_SP_8                  =   find(xwork_SP_8>0.50,1);
SP_M_8                   =   xout_SP_8(ix_SP_8);


[n_SP_BF_8,xout_SP_8_BF]    =   hist(log10(s_SP8_BF),nbin);
[xmin_SP_BF_8,tmin_SP_8_BF] =   find(n_SP_BF_8 == max(n_SP_BF_8),1);
SP_BF_8                     =   xout_SP_8_BF(tmin_SP_8_BF);

% For mu
[n_mu_8,xout_mu_8]       =   hist(s_mu8,nbin);
xwork_mu_8               =   cumsum(n_mu_8/sum(n_mu_8));

ix_mu_8                  =   find(xwork_mu_8>0.175,1);
mu_1sd_8                 =   xout_mu_8(ix_mu_8);

ix_mu_8                  =   find(xwork_mu_8>0.825,1);
mu_1su_8                 =   xout_mu_8(ix_mu_8);

ix_mu_8                  =   find(xwork_mu_8>0.025,1);
mu_2sd_8                 =   xout_mu_8(ix_mu_8);

ix_mu_8                  =   find(xwork_mu_8>0.925,1);
mu_2su_8                 =   xout_mu_8(ix_mu_8);

ix_mu_8                  =   find(xwork_mu_8>0.50,1);
mu_M_8                   =   xout_mu_8(ix_mu_8);

[n_mu_BF_8,xout_mu_8_BF]    =   hist(s_mu8_BF,nbin);
[xmin_mu_8_BF,tmin_mu_8_BF] =   find(n_mu_BF_8 == max(n_mu_BF_8),1);
mu_BF_8                     =   xout_mu_8_BF(tmin_mu_8_BF);

%% Print the results

fprintf('\nResult for the calibration with only MBTP7 alone \n')

SP_M_7_ok   = 10^(SP_M_7);
SP_BF_7_ok  = 10^(SP_BF_7);
SP_1su_7_ok = 10^(SP_1su_7);
SP_1sd_7_ok = 10^(SP_1sd_7);
SP_2su_7_ok = 10^(SP_2su_7);
SP_2sd_7_ok = 10^(SP_2sd_7);

disp(['SP Median     = ' num2str(SP_M_7_ok,3) ' s-1']);
disp(['SP BestFit    = ' num2str(SP_BF_7_ok,3) ' s-1']);
disp(['SP 1sigma sup = ' num2str(SP_1su_7_ok,3) ' s-1']);
disp(['SP 1sigma inf = ' num2str(SP_1sd_7_ok,3) ' s-1']);
disp(['SP 2sigma sup = ' num2str(SP_2su_7_ok,3) ' s-1']);
disp(['SP 2sigma inf = ' num2str(SP_2sd_7_ok,3) ' s-1']);

disp(['mu Median     = ' num2str(mu_M_7,3) ' mm-1']);
disp(['mu BestFit    = ' num2str(mu_BF_7,3) ' mm-1']);
disp(['mu 1sigma sup = ' num2str(mu_1su_7,3) ' mm-1']);
disp(['mu 1sigma inf = ' num2str(mu_1sd_7,3) ' mm-1']);
disp(['mu 2sigma sup = ' num2str(mu_2su_7,3) ' mm-1']);
disp(['mu 2sigma inf = ' num2str(mu_2sd_7,3) ' mm-1']);


fprintf('\nResult for the calibration with only MBTP8 alone \n')

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

fprintf('\n Result for the combined calibration with MBTP7 and MBTP8 \n')

SP_M_all_ok   = 10^(SP_M_all);
SP_BF_all_ok  = 10^(SP_BF_all);
SP_1su_all_ok = 10^(SP_1su_all);
SP_1sd_all_ok = 10^(SP_1sd_all);
SP_2su_all_ok = 10^(SP_2su_all);
SP_2sd_all_ok = 10^(SP_2sd_all);

disp(['SP Median     = ' num2str(SP_M_all_ok,3) ' s-1']);
disp(['SP BestFit    = ' num2str(SP_BF_all_ok,3) ' s-1']);
disp(['SP 1sigma sup = ' num2str(SP_1su_all_ok,3) ' s-1']);
disp(['SP 1sigma inf = ' num2str(SP_1sd_all_ok,3) ' s-1']);
disp(['SP 2sigma sup = ' num2str(SP_2su_all_ok,3) ' s-1']);
disp(['SP 2sigma inf = ' num2str(SP_2sd_all_ok,3) ' s-1']);

disp(['mu Median     = ' num2str(mu_M_all,3) ' mm-1']);
disp(['mu BestFit    = ' num2str(mu_BF_all,3) ' mm-1']);
disp(['mu 1sigma sup = ' num2str(mu_1su_all,3) ' mm-1']);
disp(['mu 1sigma inf = ' num2str(mu_1sd_all,3) ' mm-1']);
disp(['mu 2sigma sup = ' num2str(mu_2su_all,3) ' mm-1']);
disp(['mu 2sigma inf = ' num2str(mu_2sd_all,3) ' mm-1']);


%% Creating output models

IR50_7_M   = exp(-SP_M_7_ok*t7*exp(-mu_M_7*x_model));
IR50_7_BF  = exp(-SP_BF_7_ok*t7*exp(-mu_BF_7*x_model));

IR50_8_M   = exp(-SP_M_8_ok*t8*exp(-mu_M_8*x_model));
IR50_8_BF  = exp(-SP_BF_8_ok*t8*exp(-mu_BF_8*x_model));

IR50_7_all_M   = exp(-SP_M_all_ok*t7*exp(-mu_M_all*x_model));
IR50_7_all_BF  = exp(-SP_BF_all_ok*t7*exp(-mu_BF_all*x_model));

IR50_8_all_M   = exp(-SP_M_all_ok*t8*exp(-mu_M_all*x_model));
IR50_8_all_BF  = exp(-SP_BF_all_ok*t8*exp(-mu_BF_all*x_model));

%% Plotting

figure('NumberTitle','off','Position',[00 00 1000 1000])
sl =5;

subplot(3,2,1)
errorbar(x7(1:n7)     ,IR50_7(1:n7)     ,IR50_e7(1:n7)     ,'go' ,'markersize', sl)
% errorbar(x7(1:n7C1)     ,IR50_7(1:n7C1)     ,IR50_e7(1:n7C1)     ,'bo' ,'markersize', sl)
hold on
% errorbar(x7(n7C1+1:n7C2),IR50_7(n7C1+1:n7C2),IR50_e7(n7C1+1:n7C2),'b>' ,'markersize', sl)
% errorbar(x7(n7C2+1:n7)  ,IR50_7(n7C2+1:n7)  ,IR50_e7(n7C2+1:n7)  ,'bd' ,'markersize', sl)

plot(x_model,IR50_7_M ,'k--', 'LineWidth',1)
plot(x_model,IR50_7_BF,'k', 'LineWidth',0.5)
plot(x_model,IR50_7_all_M ,'r--', 'LineWidth',1)
plot(x_model,IR50_7_all_BF,'r', 'LineWidth',0.5)

xlabel('Depth [mm]')
ylabel('IRSL Normalized Intensity')
% legend('Core 1','Core 2','Core 3','MBTP7 Median','MBTP7 BestFit','Combined Median','Combined BestFit','Location','Southeast')
legend('Experiemntal data','MBTP7 Median','MBTP7 BestFit','Combined Median','Combined BestFit','Location','Southeast')
title('MBTP7 - IR50')
axis([0 28 0 1.3])

subplot(3,2,3)
errorbar(x8(1:n8)     ,IR50_8(1:n8)     ,IR50_e8(1:n8)     ,'go' ,'markersize', sl); 
% errorbar(x8(1:n8C1)     ,IR50_8(1:n8C1)     ,IR50_e8(1:n8C1)     ,'go' ,'markersize', sl)
 hold on
% errorbar(x8(n8C1+1:n8C2),IR50_8(n8C1+1:n8C2),IR50_e8(n8C1+1:n8C2),'g>' ,'markersize', sl)
% errorbar(x8(n8C2+1:n8)  ,IR50_8(n8C2+1:n8)  ,IR50_e8(n8C2+1:n8)  ,'gd' ,'markersize', sl)

plot(x_model,IR50_8_M ,'k--', 'LineWidth',1)
plot(x_model,IR50_8_BF,'k', 'LineWidth',0.5)
plot(x_model,IR50_8_all_M ,'r--', 'LineWidth',1)
plot(x_model,IR50_8_all_BF,'r', 'LineWidth',0.5)

xlabel('Depth [mm]')
ylabel('IRSL Normalized Intensity')
%legend('Core 1','Core 2','Core 3','MBTP8 Median','MBTP8 BestFit','Combined Median','Combined BestFit','Location','Southeast')
legend('Experiemntal data','MBTP8 Median','MBTP8 BestFit','Combined Median','Combined BestFit','Location','Southeast')
title('MBTP8 - IR50')
axis([0 28 0 1.3])


SP_vec_plot         = 0:0.01:4;
SP_M_8_line_plot    = SP_M_8_ok.*ones(size(SP_vec_plot));
SP_BF_8_line_plot   = SP_BF_8_ok.*ones(size(SP_vec_plot));
SP_M_7_line_plot    = SP_M_7_ok.*ones(size(SP_vec_plot));
SP_BF_7_line_plot   = SP_BF_7_ok.*ones(size(SP_vec_plot));
SP_M_all_line_plot  = SP_M_all_ok.*ones(size(SP_vec_plot));
SP_BF_all_line_plot = SP_BF_all_ok.*ones(size(SP_vec_plot));

mu_vec_plot         = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3];
mu_M_8_line_plot    = mu_M_8.*ones(size(mu_vec_plot));
mu_BF_8_line_plot   = mu_BF_8.*ones(size(mu_vec_plot));
mu_M_7_line_plot    = mu_M_7.*ones(size(mu_vec_plot));
mu_BF_7_line_plot   = mu_BF_7.*ones(size(mu_vec_plot));
mu_M_all_line_plot  = mu_M_all.*ones(size(mu_vec_plot));
mu_BF_all_line_plot = mu_BF_all.*ones(size(mu_vec_plot));

subplot(3,2,2)
surface((r_SP),r_mu,norm_chi7); axis square; colorbar;shading interp;set(gca,'XScale','log');
shading interp
hold on
plot(mu_vec_plot,mu_M_7_line_plot,'w--', 'LineWidth',0.5)
plot(SP_M_7_line_plot,SP_vec_plot,'w--', 'LineWidth',0.5)
plot(mu_vec_plot,mu_BF_7_line_plot,'w', 'LineWidth',0.2)
plot(SP_BF_7_line_plot,SP_vec_plot,'w', 'LineWidth',0.2)


xlabel('Bleaching rate [s-1]')
ylabel('Attenuation coeff. [mm-1]')
axis([1e-10 1e-3 0.1 4])
title('Cal. MBTP7 - IR50')

subplot(3,2,4)
surface(r_SP,r_mu,norm_chi8); axis square; colorbar;shading interp;set(gca,'XScale','log');
hold on

plot(mu_vec_plot,mu_M_8_line_plot,'w--', 'LineWidth',0.5)
plot(SP_M_8_line_plot,SP_vec_plot,'w--', 'LineWidth',0.5)
plot(mu_vec_plot,mu_BF_8_line_plot,'w', 'LineWidth',0.2)
plot(SP_BF_8_line_plot,SP_vec_plot,'w', 'LineWidth',0.2)


xlabel('Bleaching rate [s-1]')
ylabel('Attenuation coeff. [mm-1]')
axis([1e-10 1e-3 0.1 4])
title('Cal. MBTP8- IR50')

subplot(3,2,6)
surface(r_SP,r_mu,norm_chi_all); axis square; colorbar;shading interp;set(gca,'XScale','log');
hold on

plot(mu_vec_plot,mu_M_all_line_plot,'w--', 'LineWidth',0.5)
plot(SP_M_all_line_plot,SP_vec_plot,'w--', 'LineWidth',0.5)
plot(mu_vec_plot,mu_BF_all_line_plot,'w', 'LineWidth',0.2)
plot(SP_BF_all_line_plot,SP_vec_plot,'w', 'LineWidth',0.2)

xlabel('Bleaching rate [s-1]')
ylabel('Attenuation coeff. [mm-1]')
axis([1e-10 1e-3 0.1 4])
title('Cal. MBTP7 and MBTP8 combined - IR50')



