function [concentration]=cosmo_conc2age(CosmoAge,P0,er_max,ts,Lambda10Be,VecMu)

% Cosmogenic nuclide dating age with time-dependant erosion
% dN/dt   = - N*lambda+P                     Equation (3) Lal, 1991
% P(z,t)  = P(0,t)*exp(-mu*(z(t)-epsilon))   Equation (4) Lal, 1991
% lambda  = desintegration constant [yr-1]
% mu      = depth attenuation [m-1]
% epsilon = erosion rate [m/yr] 
% z       = depth [m]

% Constant of comsogenic nuclide dating
% lambda  = 4.997456e-07; % [yr-1] 
% mu      = 2.7/160*100 ; % [m-1]

% Setting of the depth referential 
zend    = 5;             %[m]
dz      = 0.01;          %[m]
z       = 0:dz:zend;     %[m]

% Numerical discretization
N2         = zeros(size(z));
dt         = min(dz/er_max,1.)/4.1;
nt         = floor(CosmoAge/dt)+1;
time       = zeros(nt,1);
Nt         = zeros(size(time));
e_rate_mat = zeros(size(time));

for it=2:nt
     time(it)       = time(it-1)+dt;
     if time(it)   <= CosmoAge-ts 
        e_rate      = 0; 
     else
        e_rate      = er_max;
     end
     advection      = e_rate*(-N2(3:end)+4*N2(2:end-1)-3*N2(1:end-2))/(2*dz);
     N2(1:end-2)    = N2(1:end-2)-dt*N2(1:end-2)*Lambda10Be+dt*P0*exp(-VecMu*z(1:end-2))+dt*advection;
     N2(N2<0)       = 0.;
     N2(end-1:end)  = N2(end-2); 
     Nt(it)         = N2(1);
     e_rate_mat(it) = e_rate;

concentration=N2(1);
% disp(P0)
% disp(Lambda10Be)
% disp(VecMu)
end
