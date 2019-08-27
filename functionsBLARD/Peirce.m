function [or_data, err_or_data] = Peirce(X,errX)
% Apply_peirce implements outlier rejection using Peirce's method as laid out by
% Stephen M. Ross, Journal of Engineering Technology, Fall 2003.  Original
% method is designed for normally distributed observations.
%
% Modified by PH Blard in May 2016 to be implemented in the CREp program
% crep.crpg.cnrs-nancy.fr
% Now include the error-weighted standard deviation mean instead of the
%
% Input : X : measured data. Must be 2 columns, n lines of measurement (column 1) with associated
% individual uncertainties at 1 sigma (column 2)
%
% Output:
% WM: New weighted mean
% WSDEV: New weighted standard deviation
% or_data: cleaned data (data without the rejected outliers)
% err_or_data : associated uncertainty
% outliers: consists of the rejected samples
% info: an optional struct that contains some information about how the process was
% carried out, useful only for the curious
%
% Example: [wheat_data chaff_data] = apply_peirce(mysampledata)
%
% Words of caution: do not use outlier rejection algorithms as a crutch;
% try to get a 'feel' for your data and use methods like these as a
% guidepost, not a rule.  Best of luck.

or_data = [];
err_or_data = [];
outliers = [];
err_outliers = [];

od = X;


% Find input dimension
%[nl,nc]=size(Data);

%
% if nl==2 && nc==2;
%     fprintf('Warning, 2 x 2 matrix, column data is the default mode\n')
% end
%
% if nc==2;
%     errX=Data(:,2);
%     X=Data(:,1);

% elseif nl==2;
%     errX=Data(2,:);
%     X=X(1,:);
% else
%     Wm=NaN;
%     fprintf('No dimension of X equal to 2\n')
%     return
% end

N=length(X);


% Load table;
peirce_r = [-1,1,2,3,4,5,6,7,8,9;3,1.19600000000000,-1,-1,-1,-1,-1,-1,-1,-1;4,1.38300000000000,1.07800000000000,-1,-1,-1,-1,-1,-1,-1;5,1.50900000000000,1.20000000000000,-1,-1,-1,-1,-1,-1,-1;6,1.61000000000000,1.29900000000000,1.09900000000000,-1,-1,-1,-1,-1,-1;7,1.69300000000000,1.38200000000000,1.18700000000000,1.02200000000000,-1,-1,-1,-1,-1;8,1.76300000000000,1.45300000000000,1.26100000000000,1.10900000000000,-1,-1,-1,-1,-1;9,1.82400000000000,1.51500000000000,1.32400000000000,1.17800000000000,1.04500000000000,-1,-1,-1,-1;10,1.87800000000000,1.57000000000000,1.38000000000000,1.23700000000000,1.11400000000000,-1,-1,-1,-1;11,1.92500000000000,1.61900000000000,1.43000000000000,1.28900000000000,1.17200000000000,1.05900000000000,-1,-1,-1;12,1.96900000000000,1.66300000000000,1.47500000000000,1.33600000000000,1.22100000000000,1.11800000000000,1.00900000000000,-1,-1;13,2.00700000000000,1.70400000000000,1.51600000000000,1.37900000000000,1.26600000000000,1.16700000000000,1.07000000000000,-1,-1;14,2.04300000000000,1.74100000000000,1.55400000000000,1.41700000000000,1.30700000000000,1.21000000000000,1.12000000000000,1.02600000000000,-1;15,2.07600000000000,1.77500000000000,1.58900000000000,1.45300000000000,1.34400000000000,1.24900000000000,1.16400000000000,1.07800000000000,-1;16,2.10600000000000,1.80700000000000,1.62200000000000,1.48600000000000,1.37800000000000,1.28500000000000,1.20200000000000,1.12200000000000,1.03900000000000;17,2.13400000000000,1.83600000000000,1.65200000000000,1.51700000000000,1.40900000000000,1.31800000000000,1.23700000000000,1.16100000000000,1.08400000000000;18,2.16100000000000,1.86400000000000,1.68000000000000,1.54600000000000,1.43800000000000,1.34800000000000,1.26800000000000,1.19500000000000,1.12300000000000;19,2.18500000000000,1.89000000000000,1.70700000000000,1.57300000000000,1.46600000000000,1.37700000000000,1.29800000000000,1.22600000000000,1.15800000000000;20,2.20900000000000,1.91400000000000,1.73200000000000,1.59900000000000,1.49200000000000,1.40400000000000,1.32600000000000,1.25500000000000,1.19000000000000;21,2.23000000000000,1.93800000000000,1.75600000000000,1.62300000000000,1.51700000000000,1.42900000000000,1.35200000000000,1.28200000000000,1.21800000000000;22,2.25100000000000,1.96000000000000,1.77900000000000,1.64600000000000,1.54000000000000,1.45200000000000,1.37600000000000,1.30800000000000,1.24500000000000;23,2.27100000000000,1.98100000000000,1.80000000000000,1.66800000000000,1.56300000000000,1.47500000000000,1.39900000000000,1.33200000000000,1.27000000000000;24,2.29000000000000,2,1.82100000000000,1.68900000000000,1.58400000000000,1.49700000000000,1.42100000000000,1.35400000000000,1.29300000000000;25,2.30700000000000,2.01900000000000,1.84000000000000,1.70900000000000,1.60400000000000,1.51700000000000,1.44200000000000,1.37500000000000,1.31500000000000;26,2.32400000000000,2.03700000000000,1.85900000000000,1.72800000000000,1.62400000000000,1.53700000000000,1.46200000000000,1.39600000000000,1.33600000000000;27,2.34100000000000,2.05500000000000,1.87700000000000,1.74600000000000,1.64200000000000,1.55600000000000,1.48100000000000,1.41500000000000,1.35600000000000;28,2.35600000000000,2.07100000000000,1.89400000000000,1.76400000000000,1.66000000000000,1.57400000000000,1.50000000000000,1.43400000000000,1.37500000000000;29,2.37100000000000,2.08800000000000,1.91100000000000,1.78100000000000,1.67700000000000,1.59100000000000,1.51700000000000,1.45200000000000,1.39300000000000;30,2.38500000000000,2.10300000000000,1.92700000000000,1.79700000000000,1.69400000000000,1.60800000000000,1.53400000000000,1.46900000000000,1.41100000000000;31,2.39900000000000,2.11800000000000,1.94200000000000,1.81200000000000,1.71000000000000,1.62400000000000,1.55000000000000,1.48600000000000,1.42800000000000;32,2.41200000000000,2.13200000000000,1.95700000000000,1.82800000000000,1.72500000000000,1.64000000000000,1.56700000000000,1.50200000000000,1.44400000000000;33,2.42500000000000,2.14600000000000,1.97100000000000,1.84200000000000,1.74000000000000,1.65500000000000,1.58200000000000,1.51700000000000,1.45900000000000;34,2.43800000000000,2.15900000000000,1.98500000000000,1.85600000000000,1.75400000000000,1.66900000000000,1.59700000000000,1.53200000000000,1.47500000000000;35,2.45000000000000,2.17200000000000,1.99800000000000,1.87000000000000,1.76800000000000,1.68300000000000,1.61100000000000,1.54700000000000,1.48900000000000;36,2.46100000000000,2.18400000000000,2.01100000000000,1.88300000000000,1.78200000000000,1.69700000000000,1.62400000000000,1.56100000000000,1.50400000000000;37,2.47200000000000,2.19600000000000,2.02400000000000,1.89600000000000,1.79500000000000,1.71100000000000,1.63800000000000,1.57400000000000,1.51700000000000;38,2.48300000000000,2.20800000000000,2.03600000000000,1.90900000000000,1.80700000000000,1.72300000000000,1.65100000000000,1.58700000000000,1.53100000000000;39,2.49400000000000,2.21900000000000,2.04700000000000,1.92100000000000,1.82000000000000,1.73600000000000,1.66400000000000,1.60000000000000,1.54400000000000;40,2.50400000000000,2.23000000000000,2.05900000000000,1.93200000000000,1.83200000000000,1.74800000000000,1.67600000000000,1.61300000000000,1.55600000000000;41,2.51400000000000,2.24100000000000,2.07000000000000,1.94400000000000,1.84300000000000,1.76000000000000,1.68800000000000,1.62500000000000,1.56800000000000;42,2.52400000000000,2.25100000000000,2.08100000000000,1.95500000000000,1.85500000000000,1.77100000000000,1.69900000000000,1.63600000000000,1.58000000000000;43,2.53300000000000,2.26100000000000,2.09200000000000,1.96600000000000,1.86600000000000,1.78300000000000,1.71100000000000,1.64800000000000,1.59200000000000;44,2.54200000000000,2.27100000000000,2.10200000000000,1.97600000000000,1.87600000000000,1.79400000000000,1.72200000000000,1.65900000000000,1.60300000000000;45,2.55100000000000,2.28100000000000,2.11200000000000,1.98700000000000,1.88700000000000,1.80400000000000,1.73300000000000,1.67000000000000,1.61400000000000;46,2.56000000000000,2.29000000000000,2.12200000000000,1.99700000000000,1.89700000000000,1.81500000000000,1.74300000000000,1.68100000000000,1.62500000000000;47,2.56800000000000,2.29900000000000,2.13100000000000,2.00600000000000,1.90700000000000,1.82500000000000,1.75400000000000,1.69100000000000,1.63600000000000;48,2.57700000000000,2.30800000000000,2.14000000000000,2.01600000000000,1.91700000000000,1.83500000000000,1.76400000000000,1.70100000000000,1.64600000000000;49,2.58500000000000,2.31700000000000,2.14900000000000,2.02600000000000,1.92700000000000,1.84400000000000,1.77300000000000,1.71100000000000,1.65600000000000;50,2.59200000000000,2.32600000000000,2.15800000000000,2.03500000000000,1.93600000000000,1.85400000000000,1.78300000000000,1.72100000000000,1.66600000000000;51,2.60000000000000,2.33400000000000,2.16700000000000,2.04400000000000,1.94500000000000,1.86300000000000,1.79200000000000,1.73000000000000,1.67500000000000;52,2.60800000000000,2.34200000000000,2.17500000000000,2.05200000000000,1.95400000000000,1.87200000000000,1.80200000000000,1.74000000000000,1.68500000000000;53,2.61500000000000,2.35000000000000,2.18400000000000,2.06100000000000,1.96300000000000,1.88100000000000,1.81100000000000,1.74900000000000,1.69400000000000;54,2.62200000000000,2.35800000000000,2.19200000000000,2.06900000000000,1.97200000000000,1.89000000000000,1.82000000000000,1.75800000000000,1.70300000000000;55,2.62900000000000,2.36500000000000,2.20000000000000,2.07700000000000,1.98000000000000,1.89800000000000,1.82800000000000,1.76700000000000,1.71100000000000;56,2.63600000000000,2.37300000000000,2.20700000000000,2.08500000000000,1.98800000000000,1.90700000000000,1.83700000000000,1.77500000000000,1.72000000000000;57,2.64300000000000,2.38000000000000,2.21500000000000,2.09300000000000,1.99600000000000,1.91500000000000,1.84500000000000,1.78400000000000,1.72900000000000;58,2.65000000000000,2.38700000000000,2.22300000000000,2.10100000000000,2.00400000000000,1.92300000000000,1.85300000000000,1.79200000000000,1.73700000000000;59,2.65600000000000,2.39400000000000,2.23000000000000,2.10900000000000,2.01200000000000,1.93100000000000,1.86100000000000,1.80000000000000,1.74500000000000;60,2.66300000000000,2.40100000000000,2.23700000000000,2.11600000000000,2.01900000000000,1.93900000000000,1.86900000000000,1.80800000000000,1.75300000000000;];

% Calc mean, std of original data

odm = mean(od);
odstd = std(od);



% Residues
orig_devs = abs(od - odm);

% populate info struct
info.odm = odm;
info.odstd = odstd;
info.orig_devs = orig_devs;
info.max_devs = [];

n = length(od);
%num_doubtful always starts at 1
num_doubtful = 1;

% find row index to use in table for this sample
n_ind = [];
if n < 61 && n > 2
    n_ind = find(peirce_r(:,1) == n);
else
    if n >= 61
        n = 60;
        warning('n > 60; using r values for n==60');
        n_ind = find(peirce_r(:,1) == n); % this line is missing from downloaded version - schauer added 140915
    end
    if n < 2
        error('this function does not support samples this small... please rethink what you are doing');
        return;
    end
end

n_ind = [];
if n < 61 && n > 2
    n_ind = find(peirce_r(:,1) == n);
else
    if n >= 61
        n = 60;
        warning('n > 60; using r values for n==60');
    end
    if n < 2
        error('this function does not support samples this small... please rethink what you are doing');
        return;
    end
end

continue_status = 1;

while(continue_status)
    
    % table_lookup
    if(num_doubtful > 9)
        error('too many outliers ( > 9) for supporting table... bailing... ');
    end
    
    r_curr = peirce_r(n_ind, num_doubtful+1);
    
    %for dbg
    %display(r_curr)
    
    % make sure valid r values in table...
    if r_curr == -1
        warning('found as many outliers as this method is capable of for this small sample... (reached end of table)');
        
        % clean up, apply outlier_rejection to data, etc
        or_data = od;
        err_or_data = errX;
        
        [nl,nc]=size(or_data);
        if nl==1;
            or_data=or_data';
        end
        
        [nl,nc]=size(err_or_data);
        if nl==1;
            err_or_data=err_or_data';
        end
        
        return;
    end
    
    %% heart of method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxdev_allowed = r_curr * odstd;
    devs = abs(od - odm);
    outlier_test = devs > maxdev_allowed;
    info.max_devs = [info.max_devs maxdev_allowed];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if tests is all zeros, then we're done,
    % if not, then we've got to go back and test for additional outliers
    
    num_outliers_found = sum(outlier_test);
    
    if (num_outliers_found > 0)
        % find outliers, reject them, and continue
        outliers = [outliers; od(outlier_test)];
        err_outliers = [err_outliers; errX(outlier_test)];
        
        % eliminate outliers
        od(outlier_test) = [];
        errX(outlier_test) = [];
        
        num_doubtful = num_doubtful + num_outliers_found;
    else
        continue_status = 0;
    end
end

% clean up
or_data = od;
err_or_data = errX;

[nl,nc]=size(or_data);
if nl==1;
    or_data=or_data';
end

[nl,nc]=size(err_or_data);
if nl==1;
    err_or_data=err_or_data';
end

return
