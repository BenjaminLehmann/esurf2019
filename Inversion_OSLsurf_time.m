%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Abs_t_one_by_one.m
% Version     : 1.02.2017 
% Specificity : division with standard deviation on the plateau of the experimental data
% Threshold   : selectionof likehood up to 0.95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc;
close all;

TT  = 10000;

% Loading data file excel MBMV_.xlsx 
Stimulation    = 'IR50_';
SampleName     = 'MBAM3';
OSL_SampleName = ['OSL_' Stimulation SampleName];
[num]          = xlsread(OSL_SampleName);
n              =  length(num(:,1));

% Input parameters from calibration
SP0 = 4.1E-06;                      % [s-1]

mu  = 0.596;                        % [mm-1]
SP  = SP0*365.25*24.*3600;          % [a-1]

%% Parameterization
% Entrer de la variable temps

tmin = 0;
tmax = 2000;   % 2*10Be age à définir

% Definitions de donnees experimentales

x(1:n)          = (num(:,1));
L(1:n)          = (num(:,2));
e(1:n)          = (num(:,3));
[x_s,ind]       = sort(x(:));
Ls              = L(ind);
a               = std(Ls(18:n));
Ddot_input      = num(1,5);
Ddot            = Ddot_input/(1e3);       % [Gy ka-1] ==> [Gy a-1]
D0              = 500;                    % [Gy]  
magicN          = Ddot/D0;                % [a-1]

%% Compute residuals (i.e. fit to data)

M         = nan(TT,1);
t_vec     = nan(TT,1);

h         = waitbar(0,'Less than one song tato...');
rand_vec  = rand(TT,1);

r_t1      = sort(tmin+(tmax-tmin)*rand_vec);


    for i = 1:TT
            
            M(i)   = 0;           
            L_th   = exp(-SP*r_t1(i)*exp(-mu*x));        % Equation without taking the dose rate in account
%            L_th   = (SP.*exp(-mu.*x).*exp(-r_t1(i).*(SP.*exp(-mu.*x)+magicN))+magicN)./(SP.*exp(-mu.*x)+magicN);    % Sohbati et al. 2012a
          %  M(i)   = sum((abs(L-L_th)./a));              % L1 norm weighted with the noise a to calculate the misfit M
            M(i)   = nansum((L-L_th).^2/a^2);             % L2 norm weighted with the noise a to calculate the misfit M         
            t_vec(i)  = r_t1(i);

        waitbar(i/TT,h)
    end
    
close(h)

chi      = 1./exp(0.5*M);    % Likelihood non normalized
max_chi  = max(chi(:));      % Max value of the Likelihood
norm_chi = chi/max_chi;      % Likelihood normalized

%%
thr      =  0.01;                        % finally we don't compare it with a random values but we take the 95% biggest
s_chi    = chi(norm_chi>thr);
s_t      = t_vec(norm_chi>thr);

%% extract 1d PDFs and confidence intervals

nbin             = 20;

[n,xout]         =   hist(s_t,nbin);
xwork            =   cumsum(n/sum(n));

ix               =   find(xwork>0.175,1);
t_1sd            =   xout(ix);
t_1sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.025,1);
t_2sd            =   xout(ix);
t_2sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.825,1);
t1_1su           =   xout(ix);
t1_su_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.925,1);
t_2su            =   xout(ix);
t_2su_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.50,1);
t_median         =   xout(ix);
t_median_pd      =   n(ix)/sum(n);

[xmin,tmin]      =   find(n == max(n),1);
t_BF             =   xout(tmin);

disp(SampleName)
T_Median  = t_median;
T_BF      = t_BF;
T_1su     = t1_1su;
T_1sd     = t_1sd;
T_2su     = t_2su;
T_2sd     = t_2sd;
%%
disp(T_Median)
%disp(T_BF)
% disp(abs(T_Median-T_1su))
% disp(abs(T_Median-T_1sd))
% disp(abs(T_Median-T_2su))
% disp(abs(T_Median-T_2sd))
disp((abs(T_Median-T_1su)+abs(T_Median-T_1sd))/2)
%save([filename '_t1.mat'],'T1_Median','T1_BF','T1_1su','T1_1sd','T1_2su','T1_2sd');
%%
xs    = 0:0.5:40;
%Ls = exp(-SP*t_median*exp(-mu*xs));
Ls_M  = (SP.*exp(-mu.*xs).*exp(-t_median.*(SP.*exp(-mu.*xs)+magicN))+magicN)./(SP.*exp(-mu.*xs)+magicN);    % Sohbati et al. 2012a
Ls_BF = (SP.*exp(-mu.*xs).*exp(-t_BF.*(SP.*exp(-mu.*xs)+magicN))+magicN)./(SP.*exp(-mu.*xs)+magicN);  

figure(102)
set(gcf,'units','points','position',[10,1200,1000,300])

subplot(1,2,1)
plot(x,L,'go','MarkerFaceColor','g')
hold on 
plot(xs,Ls_M,'r','LineWidth',1)
plot(xs,Ls_BF,'k--','LineWidth',1)
xlabel('Depth [mm]')
ylabel('Normalized IRSL Signal')
legend('Experimental values','Inversed solution Median','Inversed solution Bestfit','Location','Southeast')
axis([0 40 0 1.2])
title(['(a)  Evolution of the IRSL50 signal for ' SampleName])

t_median_vec = t_median*ones(100,1);
t_BF_vec     = t_BF*ones(100,1);
likeH        = 0:1/(100-1):1;

subplot(1,2,2)
set(gca,'XScale','log')   
plot(t_vec,norm_chi,'b','LineWidth',1)
hold on
plot(t_median_vec,likeH,'r','LineWidth',1)
plot(t_BF_vec,likeH,'k--','LineWidth',1)
axis([0 exp(tmax) 0 1])
ylabel('Likelihood')
xlabel('Time [a]')
title(['(b)  OSL surface exposur dating inversion results for Sample ' SampleName])
legend('Likelihood distribution','Median','Location','Northeast')

str1= {'OSL-surf age',...
      ['t = ' num2str(t_median,3) ' ± ' num2str(abs(T_Median-T_1sd),3) ' a']}; 
xt1 = [27 27];
yt1 = [0.475 0.4];
text(xt1,yt1,str1)
