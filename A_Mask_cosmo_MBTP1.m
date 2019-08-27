%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : A_mask_cosmo_MBTP1.m
%
% Version     : 17.07.2019
%
% Coder       : Benjamin Lehmann (lehmann.benj@gmail.com)
%
% Aim         : From 10Be concentration and considering different erosion
%               histories evolving as step function (ts: erosion time; e_r: erosion rate)
%               this code creates different matrix of size TTxTT (for each ts, e_r pairs)
%               with corrected exposure time, exposure time errors and convergence
%
% Input       : TT = matrix size
%               Sample name and 10Be dataset
%               10^log_e_max and 10^log_e_max   = Minimum and Maximum erosion rate to explore [m/a]              
%               10^log_ts_max and 10^log_ts_max = Minimum and Maximum erosion time to explore [a]
%               
% Output      : e_matrix    = Collects all the tried erosion rate 
%               ts_matrix   = Collects all the tried erosion time ts 
%               Mask_cosmo  = Convergence matrix, for each ts, e_r pairs 
%                             0 when it is not possible to calculate exposure age
%                             1 when it is possible to calculate exposure age
%               t_corr_mat  = Collects the exposure age for each ts, e_r pairs
%               Err_mat     = Collects the exposure age error for each ts, e_r pairs
%               Err2_mat    = Collects the exposure age error without PR error for each ts, e_r pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc;  
close all;

% Set matlab passes
% addpath('...');

tic()
SampleName     = 'MBTP1';
Be_SampleName  = ['Be_' SampleName];

%% Calling BeCREp to calculate exposure time assuming no erosion

[t0,Err0,Err20,Erosion_rate0,Convergence0]=BeCREp(Be_SampleName,0,0,0,0);

%% Discretization and creation , time steps and plotting frequency

TT          = 5;                   % Matrix size TT*TT
log_e_max   = -2;                  % Superior limit to explore erosion rate 
log_e_min   = -8;                  % Lower limit to explore erosion rate 
log_ts_max  = 4.5;                 % Superior limit of erosion time, to choose regarding the setting of the study (here not older than ~30 ka)
log_ts_min  = -1;                  % Lower limit to explore erosion time 
ts_r        = 10.^(linspace(log_ts_min,log_ts_max,TT));    % [a]
e_r         = 10.^(linspace(log_e_min,log_e_max,TT));      % [m/a]
 
Mask_cosmo  = nan(TT,TT);
e_matrix    = nan(TT,TT);
ts_matrix   = nan(TT,TT);
t_corr_mat  = nan(TT,TT);
Err_mat     = nan(TT,TT);
Err2_mat    = nan(TT,TT);

%% Calculation
% For Cluster work
% pool = parpool(20);

disp('The calcul is strating NOW')

for i = 1:TT                                         % Loop over e_r
disp(['Outer loop # ' num2str(i)])                   % To check where the calcule is
parfor j = 1:TT                                      % Loop over tau_r
   
    [t,Err,Err2,Erosion_rate,Convergence] = BeCREp(Be_SampleName,1,0,e_r(i),ts_r(j))
    
    e_matrix(i,j)     = e_r(i);
    ts_matrix(i,j)    = ts_r(j);
    t_corr_mat(i,j)   = t;
    Err_mat(i,j)      = Err;
    Err2_mat(i,j)     = Err2;
    
    if Convergence      == 1
        Mask_cosmo(i,j) = 1;    
    else 
        Mask_cosmo(i,j) = 0;
        
   end
end
end

%% Saving and finish

disp('Saving for the results')
name_mat = [SampleName '_TCN_Mask.mat'];
save(name_mat)
disp('Saving successfuly proceed');
toc()

% If pool where use --> necessity to close it
% delete(pool)
