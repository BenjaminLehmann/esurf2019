clear all;
close all;
clc; 

tic()

% Load data

filename      = 'MBTP1';
Be_SampleName = ['Be_' filename];
load([filename '_inversion_result.mat'])

%% Creation of the reference signal

xmax    = 0.04;                            % [m]
[t0,Err0,Err20,Erosion_rate0,Convergence0]=BeCREp(Be_SampleName,0,0,0,0);

n       = 101;
dx      = xmax/(n-1);
x_final = 0:dx:xmax;
x0      = x_final;
L0      = (sp.*exp(-mu.*x0).*exp(-t0.*(sp.*exp(-mu.*x0)+Ddot/D0))+Ddot/D0)./(sp.*exp(-mu.*x0)+Ddot/D0);    % Sohbati et al. 2012a 
L00     = exp(-sp*t0*exp(-mu*x0));

%% Plot of the data

figure(1)
set(gcf,'units','points','position',[10,1200,1200,400])

subplot(1, 2, 1)
%errorbar(x_data,L_data,L_er_data,'go','MarkerFaceColor','r','MarkerFaceColor','g')
plot(x_data,L_data,'go','MarkerFaceColor','r','MarkerFaceColor','g')
hold on
plot(x0,L0,'k--','LineWidth',2)
axis([0 0.04 0 1.2])

%% Selection of the maximum likelihood

[C,D]     = find(norm_chi>0.99);
for i     = 1:numel(C)
    t_int_lim = max(max(t_corr_mat(C(i),D(i))));
    Err_max   = max(max(Err_mat(C(i),D(i))));
    Err2_max  = max(max(Err2_mat(C(i),D(i))));
end
disp(numel(C))

%% Create synthetic luminescence with the best solution ts and e_r

for i=1:1
dt      = 0.05;                 % Here time step can be bigger than during calclulation
nnt     = floor(t_int_lim/dt);
xs      = x_final;
Ls      = ones(1,n); 
Ls(1,1) = 0;
time    = zeros(nnt,1);
time(1) = 0;
Lsp     = Ls;
res     = 0;
%
for it          = 2:nnt
    time(it)    = time(it-1)+dt;

    if time  <= t_corr_mat(C(i),D(i))-ts_matrix(C(i),D(i)) 
        s_e_rate   = 0;
    else
        s_e_rate   = e_matrix(C(i),D(i));
    end 
    advection      = s_e_rate*(-Ls(3:end)+4*Ls(2:end-1)-3*Ls(1:end-2))/(2*dx);
    Ls(1:n-2)      = Ls(1:n-2)+dt.*(Ddot/D0*(1-Ls(1:n-2))-(Ls(1:n-2)*sp).*exp(-mu.*xs(1:n-2)))+advection*dt;
    Ls             = max(Ls,0);
    Ls(end-1:end)  = 1;
    
    res = sum(abs(Lsp-Ls));
    if (res<1e-8)
        disp(['Steady state for solution = ' num2str(time(it)) ' a'])
    break
    end
    Lsp = Ls;
end

ind_min  = min(find(xs>0));         % search the index when x strats to be positive 
ind_max  = min(find(xs>0.04));      % search the index when x strats to be positive 
xs       = xs(1,ind_min:end); 
Ls       = Ls(1,ind_min:end);

subplot(1, 2, 1)
plot(xs,Ls,'r','LineWidth',1.5)
drawnow 
end

subplot(1, 2, 1)
title(['(a)  Evolution of the signal ' filename])
xlabel('Depth [m]')
ylabel('Normalized Luminescence Signal')
legend('Experimental values','Model w/o erosion','Inversed solution ind','Inversed solution all','Location','southeast')

%% Plot Porbability distribution

subplot(1, 2, 2)
set(gca,'YScale','log')
set(gca,'XScale','log')   
surface(ts_matrix,e_matrix,norm_chi)
hold on
shading interp 
% plot(ts_r,Erosion_rate_vec,'r','Linewidth',2)
axis square
axis([0.5 max(ts_r) min(e_r) max(e_r)])
colorbar 
box on
xlabel('ts [a]')
ylabel([ char(949) ' [m a^{-1}]'])
title(['(d)  Probability distribution' filename])
axis square

%% Saving

name_fig = [filename '_Inversion.fig'];
name_pdf = [filename '_Inversion.pdf'];
saveas(gcf,name_fig)
saveas(gcf,name_pdf)
toc()