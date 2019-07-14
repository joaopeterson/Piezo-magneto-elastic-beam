
% -----------------------------------------------------------------
%  main_piezomagbeam_poincare.m
%
%  This is the main file for a program that computes a Poincare
%  map for the nonlinear dynamics of a piezo-magneto-elastic beam,
%  which evolves according to the follwing system of
%  ordinary differential equations
%
%    d2x/dt2 + 2*ksi*dx/dt - 0.5*x*(1-x^2) - chi*v = f*cos(Omega*t)
%
%    dv/dt + lambda*v + kappa*dx/dt = 0
%
%        +
%
%    initial conditions,
%  
%  where
%  
%   x(t)   - dimensionless displacement of the beam tip
%   v(t)   - dimensionless voltage across the load resistance
%   t      - dimensionless time
%   ksi    - mechanical damping ratio
%   chi    - dimensionless piezoeletric coupling term (mechanical)
%   f      - dimensionless excitation amplitude
%   Omega  - dimensionless excitation frequency
%   lambda - dimensionless time constant reciprocal
%   kappa  - dimensionless piezoeletric coupling term (eletrical)
%  
%  Reference:
%  
%  A. Erturk, J. Hoffmann, and D. J. Inman
%  A piezomagnetoelastic structure for broadband vibration
%  energy harvesting
%  Applied Physics Letters
%  vol. 94 pp. 254102, 2009
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Dec 8, 2016
% -----------------------------------------------------------------


clc
clear all
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (Poincare map calculation)                         ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
%case_name = 'piezomagbeam_poincare_f010_w08';
case_name = 'poincare_f083_w08';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------


for i=1:2
% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
f      = 0.083;  % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency

if (i==2)
    f = 0.115;
    case_name = 'poincare_f0115_w08';
end

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% -----------------------------------------------------------



% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');


% physical parameters vector
phys_param = [ksi chi f Omega lambda kappa x0 xdot0 v0];

% initial conditions
IC = [x0; xdot0; v0];

% period of a forcing cycle
T = 2*pi/Omega;

% number of forcing cycles
Nf = 13000;

% initial dimensionless time
t0 = 0.0;

% final dimensionless time
t1 = t0 + Nf*T;

% number of samples per forcing cycle
Nsamp = 1000;

% time series sampling points
tspan = t0:(T/Nsamp):t1;

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

% ODE solver Runge-Kutta45
[time,y] = ode45(@(t,y)piezomagelastbeam_nonlinear(t,y,phys_param),tspan,IC,opt);

% time series of dimensionless displacement
Qdisp = y(:,1);

% time series of dimensionless velocity
Qvelo = y(:,2);

% time series of dimensionless voltage
Qvolt = y(:,3);

% number of dimensionless time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.99*Ndt);

% number of steps to initiates Poincare map
Npm = round(0.25*Ndt);

toc
% -----------------------------------------------------------


% compute the Poincare map
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- computing the Poincare map --- ');
disp(' ');
disp('    ... ');
disp(' ');

% Poincare maps
poincare_disp = Qdisp(Npm:Nsamp:Ndt);
poincare_velo = Qvelo(Npm:Nsamp:Ndt);
poincare_volt = Qvolt(Npm:Nsamp:Ndt);

toc
% -----------------------------------------------------------



% save simulation results
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------



% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');



% plot disp velo Poincare map
% ...........................................................
gtitle = '';
xlab   = ' displacement';
ylab   = ' velocity';
xmin   = -2;
xmax   =  2;
ymin   = -2;
ymax   =  2;
gname  = [num2str(case_name),'__disp_velo'];
flag   = 'png';
fig1   = plot_poincare_map(Qdisp(Nss:Ndt),Qvelo(Nss:Ndt),...
                           poincare_disp,poincare_velo,gtitle,....
                           xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ...........................................................



% plot displacement Poincare map
% ...........................................................
gtitle = '';
xlab   = ' displacement';
ylab   = ' voltage';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
xmin   = -2;
xmax   =  2;
ymin   = -1;
ymax   =  1;
gname  = [num2str(case_name),'__disp_volt'];
flag   = 'png';
fig2   = plot_poincare_map(Qdisp(Nss:Ndt),Qvolt(Nss:Ndt),...
                           poincare_disp,poincare_volt,gtitle,....
                           xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2);
% ...........................................................


% plot displacement Poincare map
% ...........................................................
gtitle = '';
xlab   = ' velocity';
ylab   = ' voltage';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
xmin   = -2;
xmax   =  2;
ymin   = -1;
ymax   =  1;
gname  = [num2str(case_name),'__velo_volt'];
flag   = 'png';
fig3   = plot_poincare_map(Qvelo(Nss:Ndt),Qvolt(Nss:Ndt),...
                           poincare_velo,poincare_volt,gtitle,....
                           xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3);
% ...........................................................

toc
% -----------------------------------------------------------
end
