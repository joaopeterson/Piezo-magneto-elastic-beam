%% -----------------------------------------------------------------
%  main_piezomagelastbeam_ivp.m
%
%  This is the main file for a program that simulates the nonlinear
%  dynamics of a piezo-magneto-elastic beam, which evolves according to the
%  follwing system of ordinary differential equations
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
%   x(t)   - dimensionless displacement of the beam tip v(t)   -
%   dimensionless voltage across the load resistance t      - dimensionless
%   time ksi    - mechanical damping ratio chi    - dimensionless
%   piezoeletric coupling term (mechanical) f      - dimensionless
%   excitation amplitude Omega  - dimensionless excitation frequency lambda
%   - dimensionless time constant reciprocal kappa  - dimensionless
%   piezoeletric coupling term (eletrical)
%  
%  Reference:
%  
%  A. Erturk, J. Hoffmann, and D. J. Inman A piezomagnetoelastic structure
%  for broadband vibration energy harvesting Applied Physics Letters vol.
%  94 pp. 254102, 2009
% -----------------------------------------------------------------
%  programmers: Joao Victor Ligier Lopes Peterson
%               joao.peterson@hotmail.com
%               
%               Vinicius Goncalves Lopes 
%               vinigolop@hotmail.com
%
%               Americo Barbosa da Cunha Junior 
%               americo.cunhajr@gmail.com
%
%  last update: Jan 12, 2016
% -----------------------------------------------------------------
%%

clc
clear
close all


%% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (nonlinear dynamics integration)                   ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Joao Victor Ligier Lopes Peterson                  ')
disp(' joao.peterson@hotmail.com                          ')
disp('                                                    ')
disp(' Vinicius Goncalves Lopes                           ')
disp(' vinigolop@hotmail.com                              ')
disp('                                                    ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
% -----------------------------------------------------------

%% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

%case_name = 'nonlinear';
case_name = 'linear_vs_nonlinear';

for a=1:9
    
ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
f      = 0.115; % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)

%f = (a-1)*0.032+0.019;
Omega = a/10;

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% physical parameters vector
phys_param = [ksi chi f Omega lambda kappa x0 xdot0 v0];
% -----------------------------------------------------------

%% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');

t0 = 0.0;                                      % initial dimensionless time
t1 = 5.0e3;                                    % final dimensionless time

IC = [x0; xdot0; v0];                          % initial conditions

opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9); % ODE solver optional parameters

% ODE solver Runge-Kutta45
[time_nonlinear,y_nonlinear] = ode45(@(t,y)piezomagelastbeam_nonlinear(t,y,phys_param),[t0 t1],IC,opt);
[time_linear,y_linear] = ode45(@(t,y)piezomagelastbeam_linear(t,y,phys_param),[t0 t1],IC,opt);

Qdisp_nonlinear = y_nonlinear(:,1);
Qvelo_nonlinear = y_nonlinear(:,2);
Qvolt_nonlinear = y_nonlinear(:,3);

Qdisp_linear = y_linear(:,1);
Qvelo_linear = y_linear(:,2);
Qvolt_linear = y_linear(:,3);

Ndt_nonlinear = length(time_nonlinear);    % number of time steps
Nss_nonlinear = round(0.9*Ndt_nonlinear); % steady state

Ndt_linear = length(time_linear);    % number of time steps
Nss_linear = round(0.98*Ndt_linear); % steady state

toc
% -----------------------------------------------------------

%% save simulation results
% -----------------------------------------------------------
% tic
% 
% disp(' ')
% disp(' --- saving simulation results --- ');
% disp(' ');
% disp('    ... ');
% disp(' ');
% 
% save('chaotic_attractor.txt', 'y_nonlinear', '-ascii')
% 
% toc
% -----------------------------------------------------------

%% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');

% time series display window
initial_percent = 0.89;
final_percent = 0.90;

%% plot displacement
%{
Njump = 1; % number of steps to jump

X = time_nonlinear(1:Njump:end);
Y = Qdisp_nonlinear(1:Njump:end);
X1 = time_linear(1:Njump:end);
Y1 = Qdisp_linear(1:Njump:end);

X3 = time_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));
Y3 = Qdisp_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));
X4 = time_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));
Y4 = Qdisp_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));

fig1 = figure('NumberTitle','off');
ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
plot(ax1, X, Y, 'b')
hold on
plot(ax1, X1, Y1, 'r')
rectangle('position',[X3(1)         min(min(Y3),min(Y4))-0.2 ...
                      X3(end)-X3(1) abs(min(min(Y3),min(Y4))-0.2)+max(max(Y3),max(Y4))+0.2],...
              'EdgeColor', 'yellow', 'LineWidth', 2.0);
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax1,'displacement', 'FontSize', 16, 'FontName', 'Helvetica');
hold off

ax3 = axes('Position',[0.35 0.70 0.25 0.21]);
plot(ax3, X3, Y3, 'b')
hold on
plot(ax3, X4, Y4, 'r')
xlim([X3(1)-2 X3(end)]);
ylim([min(min(Y3),min(Y4))-0.15 max(max(Y3),max(Y4))+0.15]);
set(ax3,'XColor', 'yellow');
set(ax3,'YColor', 'yellow');
set(ax3,'LineWidth', 2.0);
set(ax3,'Box','on');
set(ax3,'XTickLabel','');
set(ax3,'YTickLabel','');
hold off

%gname = ['disp_time_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
%saveas(gcf, gname, 'png')
%close(fig1);
%}

%% plot velocity
%{
Njump = 1; % number of steps to jump

X = time_nonlinear(1:Njump:end);
Y = Qvelo_nonlinear(1:Njump:end);
X2 = time_linear(1:Njump:end);
Y2 = Qvelo_linear(1:Njump:end);

%initial_percent = 0.85;
%final_percent = 0.95;
X3 = time_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));
Y3 = Qvelo_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));
X4 = time_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));
Y4 = Qvelo_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));

fig2 = figure('NumberTitle','off');
ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
plot(ax1, X, Y, 'b')
hold on
plot(ax1, X2, Y2, 'r')
rectangle('position',[X3(1)         min(min(Y3),min(Y4))-0.2 ...
                      X3(end)-X3(1) abs(min(min(Y3),min(Y4))-0.2)+max(max(Y3),max(Y4))+0.2],...
              'EdgeColor', 'yellow', 'LineWidth', 2.0);
xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax1,'velocity', 'FontSize', 16, 'FontName', 'Helvetica');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
hold off

ax3 = axes('Position',[0.35 0.70 0.25 0.21]);
plot(ax3, X3, Y3, 'b')
hold on
plot(ax3, X4, Y4, 'r')
xlim([X3(1)-2 X3(end)]);
ylim([min(min(Y3),min(Y4))-0.15 max(max(Y3),max(Y4))+0.15]);
set(ax3,'XColor', 'yellow');
set(ax3,'YColor', 'yellow');
set(ax3,'LineWidth', 2.0);
set(ax3,'Box','on');
set(ax3,'XTickLabel','');
set(ax3,'YTickLabel','');
hold off

gname = ['velo_time_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
saveas(gcf, gname, 'png')
close(fig2);
% ...........................................................
%}

%% plot voltage

Njump = 1; % number of steps to jump

X = time_nonlinear(1:Njump:end);
Y = Qvolt_nonlinear(1:Njump:end);

if(strcmp(case_name,'linear_vs_nonlinear'))
    X2 = time_linear(1:Njump:end);
    Y2 = Qvolt_linear(1:Njump:end);

    X4 = time_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));
    Y4 = Qvolt_linear(round(initial_percent*Ndt_linear) : Njump : round(final_percent*Ndt_linear));
end


%initial_percent = 0.85;
%final_percent = 0.95;
%X3 = time_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));
%Y3 = Qvolt_nonlinear(round(initial_percent*Ndt_nonlinear) : Njump : round(final_percent*Ndt_nonlinear));

X3 = time_nonlinear(round(0.95*Ndt_nonlinear-(600/Omega)) : 1 : round(0.95*Ndt_nonlinear));
Y3 = Qvolt_nonlinear(round(0.95*Ndt_nonlinear-(600/Omega)) : 1 : round(0.95*Ndt_nonlinear));

fig3 = figure('NumberTitle','off');
ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
plot(ax1, X, Y, 'b')
set(gcf,'Color','White');

if(strcmp(case_name,'linear_vs_nonlinear'))
    hold on
    plot(ax1, X2, Y2, 'r')
    rectangle('position',[X3(1)         min(min(Y3),min(Y4))-0.2 ...
                          X3(end)-X3(1) abs(min(min(Y3),min(Y4))-0.2)+max(max(Y3),max(Y4))+0.2],...
                  'EdgeColor', 'yellow', 'LineWidth', 2.0);
else
    rectangle('position',[X3(1)         min(Y3)-0.2 ...
                          X3(end)-X3(1) abs(min(Y3)-0.2)+max(Y3)+0.2],...
                  'EdgeColor', 'yellow', 'LineWidth', 2.0);
end
xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax1,'voltage', 'FontSize', 16, 'FontName', 'Helvetica');
ylim([-1.2 1.2]);
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
hold off

ax3 = axes('Position',[0.55 0.70 0.25 0.21]);
plot(ax3, X3, Y3, 'b')

if(strcmp(case_name,'linear_vs_nonlinear'))
    hold on
    plot(ax3, X4, Y4, 'r')
    
    xlim([X3(1)-2 X3(end)]);
    ylim([min(min(Y3),min(Y4))-0.15 max(max(Y3),max(Y4))+0.15]);
else
    xlim([X3(1)-2 X3(end)]);
    ylim([min(Y3)-0.15 max(Y3)+0.15]);
end


set(ax3,'XColor', 'yellow');
set(ax3,'YColor', 'yellow');
set(ax3,'LineWidth', 2.0);
set(ax3,'Box','on');
set(ax3,'XTickLabel','');
set(ax3,'YTickLabel','');
hold off

if(strcmp(case_name,'linear_vs_nonlinear'))
    if(f<0.1)
        gname     = ['volt_time_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
    else
        gname     = ['volt_time_f', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
    end
else
    if(f<0.1)
        gname     = ['volt_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname     = ['volt_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end
end

saveas(gcf, gname, 'png')
close(fig3);
%}
% ...........................................................

%% plot disp vs velo
%{
X = Qdisp_nonlinear(Nss_nonlinear:Ndt_nonlinear);
Y = Qvelo_nonlinear(Nss_nonlinear:Ndt_nonlinear);

fig4 = figure('NumberTitle','off');
ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
fh1 = plot(ax1, X, Y, 'b');
hold on

if(strcmp(case_name,'linear_vs_nonlinear'))
    X1 = Qdisp_linear(Nss_linear:Ndt_linear);
    Y1 = Qvelo_linear(Nss_linear:Ndt_linear);
    fh2 = plot(ax1, X1, Y1, 'r');
    gname = ['disp_velo_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
    
elseif(strcmp(case_name,'nonlinear'))
    gname = ['disp_velo_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];
    
end

xlim(ax1, [-2 2]);
ylim(ax1,[-1.5 1.5]);
set(fh1,'LineWidth',1.0);
set(gcf,'color','white');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'displacement', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel(ax1,'velocity', 'FontSize', 18, 'FontName', 'Helvetica');
hold off


saveas(gcf, gname, 'png')
%close(fig4);
%}
% ...........................................................

%% plot disp vs volt
%{
X = Qdisp_nonlinear(Nss_nonlinear:Ndt_nonlinear);
Y = Qvolt_nonlinear(Nss_nonlinear:Ndt_nonlinear);

fig4 = figure('NumberTitle','off');
ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
fh1 = plot(ax1, X, Y, 'b');
hold on

if(strcmp(case_name,'linear_vs_nonlinear'))
    X1 = Qdisp_linear(Nss_linear:Ndt_linear);
    Y1 = Qvolt_linear(Nss_linear:Ndt_linear);
    fh2 = plot(ax1, X1, Y1, 'r');
    gname = ['disp_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
    
elseif(strcmp(case_name,'nonlinear'))
    gname = ['disp_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];
    
end

xlim(ax1, [-2 2]);
ylim(ax1,[-1.5 1.5]);
set(fh1,'LineWidth',1.0);
set(gcf,'color','white');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'displacement', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel(ax1,'voltage', 'FontSize', 18, 'FontName', 'Helvetica');
hold off


saveas(gcf, gname, 'png')
%close(fig4);
%}
% ...........................................................

%% plot velo vs volt
%{
X = Qvelo_nonlinear(Nss_nonlinear:Ndt_nonlinear);
Y = Qvolt_nonlinear(Nss_nonlinear:Ndt_nonlinear);

fig4 = figure('NumberTitle','off');
ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
fh1 = plot(ax1, X, Y, 'b');
hold on

if(strcmp(case_name,'linear_vs_nonlinear'))
    X1 = Qvelo_linear(Nss_linear:Ndt_linear);
    Y1 = Qvolt_linear(Nss_linear:Ndt_linear);
    fh2 = plot(ax1, X1, Y1, 'r');
    gname = ['velo_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];
    
elseif(strcmp(case_name,'nonlinear'))
    gname = ['velo_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];
    
end

xlim(ax1, [-2 2]);
ylim(ax1,[-1.5 1.5]);
set(fh1,'LineWidth',1.0);
set(gcf,'color','white');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'velocity', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel(ax1,'voltage', 'FontSize', 18, 'FontName', 'Helvetica');
hold off


saveas(gcf, gname, 'png')
%close(fig4);
%}
% ...........................................................

%% plot atractor
%{
fig5 = figure('NumberTitle','off');

X1 = Qdisp_nonlinear(Nss_nonlinear:Ndt_nonlinear);
Y1 = Qvelo_nonlinear(Nss_nonlinear:Ndt_nonlinear);
Z1 = Qvolt_nonlinear(Nss_nonlinear:Ndt_nonlinear);

fh1 = plot3(X1, Y1, Z1, 'b');
hold on

if(strcmp(case_name,'linear_vs_nonlinear'))
    X2 = Qdisp_linear(Nss_linear:Ndt_linear);
    Y2 = Qvelo_linear(Nss_linear:Ndt_linear);
    Z2 = Qvolt_linear(Nss_linear:Ndt_linear);
    fh2 = plot3(X2, Y2, Z2, 'r');
    gname = ['3d_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear_vs_nonlinear'];

elseif(strcmp(case_name,'nonlinear'))
    gname = ['3d_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];

end

view(45, 45);
grid on;
xlim([-2 2])
ylim([-2 2])
zlim([-1 1])
set(fh1,'LineWidth',1.0);
set(gcf,'color','white');
set(gca,'TickDir','out','TickLength',[.02 .02]);
set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'zColor',[.3 .3 .3]);
xlabel('displacement', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel('velocity', 'FontSize', 18, 'FontName', 'Helvetica');
zlabel('voltage', 'FontSize', 18, 'FontName', 'Helvetica');
hold off

saveas(gcf, gname, 'png')
%close(fig5);
%}
% ...........................................................

%% plot histograms
%{
%Nbins = round(sqrt(Ndt_nonlinear));
Nbins = 64;
Nksd  = 32*Nbins;

Ndt = length(Qvolt_nonlinear);
Nss = round(0.9*Ndt);
Pot_nonlinear = lambda.*Qvolt_nonlinear(Nss:Ndt).*Qvolt_nonlinear(Nss:Ndt);

[Pot_bins,Pot_freq] = randvar_pdf(Pot_nonlinear,Nbins);
[Pot_ksd ,Pot_supp] = randvar_ksd(Pot_nonlinear,Nksd);

gtitle    = ' ';
xlab      = ' power';
ylab      = ' probability density';
xmin      = 0.0;
xmax      = max(Pot_bins);
ymin      = 0.0;
ymax      = 'auto';
if(f<0.1)
    gname     = ['pdf_pot_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];
else
    gname     = ['pdf_pot_f', num2str(f*1000), '_O0', num2str(Omega*10), '_nonlinear'];
end
linearity = 'nonlinear';
flag      = 'png';
fig6a     = graph_bar_curve1(Pot_bins,Pot_freq,...
                          Pot_supp,Pot_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,linearity,flag);
        
%close(fig6a);

%Nbins = round(sqrt(Ndt_linear));
Nbins = 64;
Nksd  = 32*Nbins;

Ndt = length(Qvolt_linear);
Nss = round(0.9*Ndt);
Pot_linear = lambda.*Qvolt_linear(Nss:Ndt).*Qvolt_linear(Nss:Ndt);

[Pot_bins,Pot_freq] = randvar_pdf(Pot_linear,Nbins);
[Pot_ksd ,Pot_supp] = randvar_ksd(Pot_linear,Nksd);

gtitle    = ' ';
xlab      = ' power';
ylab      = ' probability density';
xmin      = 0;
xmax      = 4e-5;
ymin      = 0;
ymax      = 'auto';
if(f<0.1)
    gname     = ['pdf_pot_f0', num2str(f*1000), '_O0', num2str(Omega*10), '_linear'];
else
    gname     = ['pdf_pot_f', num2str(f*1000), '_O0', num2str(Omega*10), '_linear'];
end
linearity = 'linear';
flag      = 'png';
fig6b     = graph_bar_curve1(Pot_bins,Pot_freq,...
                          Pot_supp,Pot_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,linearity,flag);
%close(fig6b);
%}

end
%%
toc
% -----------------------------------------------------------

