% main_z1_analysis.m
%
% plot graphs of K x f / k x omega

clc
clear all
close all

% program header
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

% simulation information
% -----------------------------------------------------------
for a=1:1
    for b=1:3
        if (a==2)
            case_name = 'K_vs_f';
        else
            case_name = 'K_vs_omega';
        end


disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------

% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

t0 = 0.0;    % initial dimensionless time
t1 = 5.0e3;  % final dimensionless time

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)

lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% initial conditions
IC = [x0; xdot0; v0];

if strcmp(case_name,'K_vs_f')

    % dimensionless excitation amplitude
    f_min = 0.019;
    f_max = 0.275;
    N_f = 500;
    f = linspace(f_min, f_max, N_f);

    % dimensionless excitation frequency
    Omega = b/10;
    
    if (Omega > 0.5)
       break 
    end

    max_i = N_f;

elseif strcmp(case_name,'K_vs_omega')

    % dimensionless excitation frequency
    omega_min = 0.1;
    omega_max = 0.9;
    N_omega = 500;
    Omega = linspace(omega_min, omega_max, N_omega);

    % dimensionless excitation amplitude
    %f = (b-1)*0.032+0.019;
    f = (b-1)*0.04+0.04;
    
    max_i = N_omega;

end
% -----------------------------------------------------------

% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

% allocate memory for K_vec
K_vec = zeros(1, max_i);

N_over = 0;

for i=1:max_i
    if strcmp(case_name,'K_vs_f')
        % physical parameters vector
        phys_param = [ksi chi f(i) Omega lambda kappa];
    elseif strcmp(case_name, 'K_vs_omega')
        % physical parameters vector
        phys_param = [ksi chi f Omega(i) lambda kappa];
    end

    % ODE solver Runge-Kutta45
    [time,yresp] = ode45(@(t,y)piezomagelastbeam_nonlinear(t,y,phys_param),[t0 t1],IC,opt);

    % time series of dimensionless voltage
    Qvolt = yresp(:,3);

    % number of dimensionless time steps
    Ndt = length(time);

    % number of steps for steady state
    Nss = Ndt-30000;

    % number of steps to jump
    Njump = 60;
    
    % 0-1 test indicator
    [K_vec(i), flag] = z1test(Qvolt(Nss:Njump:Ndt));
    
    N_over = N_over + flag;
end

msg = ['Number of warnings: ', num2str(N_over)]; 
disp(msg);
disp(' ');

msg = ['Percentage of total: ', num2str(N_over*100/max_i), '%'];
disp(msg);
disp(' ');


toc
% -----------------------------------------------------------

% save simulation results
% -----------------------------------------------------------

disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

if strcmp(case_name,'K_vs_f')
    save([num2str(case_name), '_O0', num2str(Omega*10), '.mat']);
elseif strcmp(case_name, 'K_vs_omega')
    save([num2str(case_name), '_f0', num2str(f*1000), '.mat']);
end
%}
% -----------------------------------------------------------

% plotting
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- plotting --- ');
disp(' ');
disp('    ... ');
disp(' ');

if strcmp(case_name,'K_vs_f')

    fig1 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.17 0.2 0.7 0.7]);
    plot(ax1,f, K_vec, 'r.');
    xlim([f_min f_max]);
    ylim([0 1]);
    set(gcf,'color','white');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    xlabel(ax1,'excitation amplitude', 'FontSize', 18, 'FontName', 'Helvetica');
    ylabel(ax1,'test 0-1 classifier', 'FontSize', 18, 'FontName', 'Helvetica');

    gname = [num2str(case_name), '_O0', num2str(Omega*10)];
    saveas(gcf, gname, 'png')
    %close(fig1);

elseif strcmp(case_name, 'K_vs_omega')

    fig1 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.17 0.2 0.7 0.7]);

    plot(ax1,Omega, K_vec, 'r.')
    xlim([omega_min omega_max]);
    ylim([0 1]);
    set(gcf,'color','white');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    xlabel(ax1,'excitation frequency', 'FontSize', 18, 'FontName', 'Helvetica');
    ylabel(ax1,'test 0-1 classifier', 'FontSize', 18, 'FontName', 'Helvetica');

    gname = [num2str(case_name), '_f0', num2str(f*1000)];
    saveas(gcf, gname, 'png')
    %close(fig1);
end
    end
end


