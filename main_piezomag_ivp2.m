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
%% Decide Plots (1 for true)
displacement = 1;
velocity     = 1;
voltage      = 1;
disp_velo    = 0;
disp_volt    = 0;
velo_volt    = 0;
attractor    = 0;
histogram    = 0;

%% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

for a=1:1
    
ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)

if a < 10
    Omega = 0.8;              % dimensionless excitation frequency
    f = (a-1)*0.032+0.019;    % dimensionless excitation amplitude
elseif a < 19
    Omega = 0.1 + (a-10)*0.1;  % dimensionless excitation frequency
    f = 0.083;                % dimensionless excitation amplitude
else
    Omega = 0.1 + (a-19)*0.1;  % dimensionless excitation frequency
    f = 0.115;                % dimensionless excitation amplitude
end
Omega = 0.8;
f = 0.115;

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% physical parameters vector
phys_param = [ksi chi f Omega lambda kappa];

% initial conditions
IC = [x0; xdot0; v0];
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

opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9); % ODE solver optional parameters

% ODE solver Runge-Kutta45
[time,y] = ode45(@(t,y)piezomagelastbeam_nonlinear(t,y,phys_param),[t0 t1],IC,opt);

Qdisp = y(:,1);
Qvelo = y(:,2);
Qvolt = y(:,3);

Ndt  = length(time);     % number of time steps
Nss  = round(0.98*Ndt);  % steady state
Nts1 = round(0.05*Ndt);  % Start of transient state for 0-1 test
Nts2 = round(0.15*Ndt);  % End of transient state for 0-1 test

toc
% -----------------------------------------------------------

%% save simulation results
% -----------------------------------------------------------
%{
tic

disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

filename = ['y_f0', num2str(f*1000), '_O0', num2str(Omega*10), '.mat'];
save(filename,'time','y')

toc
%}
% -----------------------------------------------------------

%% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');

%% plot displacement
if displacement == 1
    Njump = 1; % number of steps to jump

    X = time(1:Njump:end);
    Y = Qdisp(1:Njump:end);
    xlab = 'time';
    ylab = 'displacement';
    
    if(f<0.1)
        gname = ['disp_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['disp_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end
    
    plot_time_series(X, Y, xlab, ylab, gname)
    %{
    X2 = time(Nss:Njump:Ndt);
    Y2 = Qdisp(Nss:Njump:Ndt);

    fig1 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
    plot(ax1, X, Y, 'b')
    ylim([-2.5 2.5]);
    hold on
    pos1 = X(Nss);
    pos2 = min(Y(Nss:Ndt))-0.1;
    pos3 = (X(end)-X(Nss))*0.2;
    pos4 = max(Y(Nss:Ndt))-min(Y(Nss:Ndt))+0.2;
    rectangle('position',[pos1, pos2, pos3, pos4],...
                  'EdgeColor', 'red', 'LineWidth', 2.0);
              
    ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
    plot(ax2, X2, Y2, 'b')
    hold on
    xlim([pos1 pos1+pos3]);
    ylim([pos2 pos2+pos4]);
    set(ax2,'XColor', 'red');
    set(ax2,'YColor', 'red');
    set(ax2,'LineWidth', 2.0);
    set(ax2,'Box','on');
    set(ax2,'XTickLabel','');
    set(ax2,'YTickLabel','');
    hold on
    
    if K > 0.8
        X3 = time(Nts1:Njump:Nts2);
        Y3 = Qdisp(Nts1:Njump:Nts2); 
        
        pos1 = X(Nts1);
        pos2 = min(Y(Nts1:Nts2))-0.1;
        pos3 = (X(Nts2)-X(Nts1))*0.2;
        pos4 = max(Y(Nts1:Nts2))-min(Y(Nts1:Nts2))+0.2;
        rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
                      'EdgeColor', 'yellow', 'LineWidth', 2.0);

        ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
        plot(ax3, X3, Y3, 'b')
        hold on
        xlim([pos1 pos1+pos3]);
        ylim([pos2 pos2+pos4]);
        set(ax3,'XColor', 'yellow');
        set(ax3,'YColor', 'yellow');
        set(ax3,'LineWidth', 2.0);
        set(ax3,'Box','on');
        set(ax3,'XTickLabel','');
        set(ax3,'YTickLabel','');
    end
    
    xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel(ax1,'displacement', 'FontSize', 16, 'FontName', 'Helvetica');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    hold off

    if(f<0.1)
        gname = ['disp_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['disp_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end

    saveas(gcf, gname, 'png')
    close(fig1);
    %}
end
% ...........................................................


%% plot velocity
if velocity == 1
    Njump = 1; % number of steps to jump

    X = time(1:Njump:end);
    Y = Qvelo(1:Njump:end);
    xlab = 'time';
    ylab = 'velocity';
    
    if(f<0.1)
        gname = ['velo_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['velo_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end
    
    plot_time_series(X, Y, xlab, ylab, gname)
    %{
    X2 = time(Nss:Njump:Ndt);
    Y2 = Qvelo(Nss:Njump:Ndt);

    fig2 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
    plot(ax1, X, Y, 'b')
    ylim([-2.5 2.5]);
    hold on
    pos1 = X(Nss);
    pos2 = min(Y(Nss:Ndt))-0.1;
    pos3 = (X(end)-X(Nss))*0.2;
    pos4 = max(Y(Nss:Ndt))-min(Y(Nss:Ndt))+0.2;
    rectangle('position',[pos1, pos2, pos3, pos4],...
                  'EdgeColor', 'red', 'LineWidth', 2.0);
              
    ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
    plot(ax2, X2, Y2, 'b')
    hold on
    xlim([pos1 pos1+pos3]);
    ylim([pos2 pos2+pos4]);
    set(ax2,'XColor', 'red');
    set(ax2,'YColor', 'red');
    set(ax2,'LineWidth', 2.0);
    set(ax2,'Box','on');
    set(ax2,'XTickLabel','');
    set(ax2,'YTickLabel','');
    hold on
    
    if K > 0.8
        X3 = time(Nts1:Njump:Nts2);
        Y3 = Qvelo(Nts1:Njump:Nts2); 
        
        pos1 = X(Nts1);
        pos2 = min(Y(Nts1:Nts2))-0.1;
        pos3 = (X(Nts2)-X(Nts1))*0.2;
        pos4 = max(Y(Nts1:Nts2))-min(Y(Nts1:Nts2))+0.2;
        rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
                      'EdgeColor', 'yellow', 'LineWidth', 2.0);

        ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
        plot(ax3, X3, Y3, 'b')
        hold on
        xlim([pos1 pos1+pos3]);
        ylim([pos2 pos2+pos4]);
        set(ax3,'XColor', 'yellow');
        set(ax3,'YColor', 'yellow');
        set(ax3,'LineWidth', 2.0);
        set(ax3,'Box','on');
        set(ax3,'XTickLabel','');
        set(ax3,'YTickLabel','');
    end
    
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

    if(f<0.1)
        gname = ['velo_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['velo_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end

    saveas(gcf, gname, 'png')
    close(fig2);
    %}
end
% ...........................................................


%% plot voltage
if voltage == 1
    Njump = 1; % number of steps to jump

    X = time(1:Njump:end);
    Y = Qvolt(1:Njump:end);
    xlab = 'time';
    ylab = 'voltage';
    
    if(f<0.1)
        gname = ['volt_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['volt_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end
    
    plot_time_series(X, Y, xlab, ylab, gname)
    %{
    X2 = time(Nss:Njump:Ndt);
    Y2 = Qvolt(Nss:Njump:Ndt);

    fig3 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
    plot(ax1, X, Y, 'b')
    ylim([-2.5 2.5]);
    hold on
    pos1 = X(Nss);
    pos2 = min(Y(Nss:Ndt))-0.1;
    pos3 = (X(end)-X(Nss))*0.2;
    pos4 = max(Y(Nss:Ndt))-min(Y(Nss:Ndt))+0.2;
    rectangle('position',[pos1, pos2, pos3, pos4],...
                  'EdgeColor', 'red', 'LineWidth', 2.0);
              
    ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
    plot(ax2, X2, Y2, 'b')
    hold on
    xlim([pos1 pos1+pos3]);
    ylim([pos2 pos2+pos4]);
    set(ax2,'XColor', 'red');
    set(ax2,'YColor', 'red');
    set(ax2,'LineWidth', 2.0);
    set(ax2,'Box','on');
    set(ax2,'XTickLabel','');
    set(ax2,'YTickLabel','');
    hold on
    
    if K > 0.8
        X3 = time(Nts1:Njump:Nts2);
        Y3 = Qvolt(Nts1:Njump:Nts2); 
        
        pos1 = X(Nts1);
        pos2 = min(Y(Nts1:Nts2))-0.1;
        pos3 = (X(Nts2)-X(Nts1))*0.2;
        pos4 = max(Y(Nts1:Nts2))-min(Y(Nts1:Nts2))+0.2;
        rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
                      'EdgeColor', 'yellow', 'LineWidth', 2.0);

        ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
        plot(ax3, X3, Y3, 'b')
        hold on
        xlim([pos1 pos1+pos3]);
        ylim([pos2 pos2+pos4]);
        set(ax3,'XColor', 'yellow');
        set(ax3,'YColor', 'yellow');
        set(ax3,'LineWidth', 2.0);
        set(ax3,'Box','on');
        set(ax3,'XTickLabel','');
        set(ax3,'YTickLabel','');
    end
    
    xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel(ax1,'voltage', 'FontSize', 16, 'FontName', 'Helvetica');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    hold off

    if(f<0.1)
        gname = ['volt_time_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['volt_time_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end

    saveas(gcf, gname, 'png')
    close(fig3);
    %}
end
% ...........................................................

%% plot disp vs velo
if disp_velo == 1
    Nss = round(0.9*Ndt);

    X = Qdisp(Nss:Ndt);
    Y = Qvelo(Nss:Ndt);

    fig4 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
    fh1 = plot(ax1, X, Y, 'b');
    hold on

    if(f<0.1)
        gname = ['disp_velo_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['disp_velo_f', num2str(f*1000), '_O0', num2str(Omega*10)];
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
    close(fig4);
end
% ...........................................................

%% plot disp vs volt
if disp_volt == 1
    Nss = round(0.9*Ndt);

    X = Qdisp(Nss:Ndt);
    Y = Qvolt(Nss:Ndt);

    fig5 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
    fh1 = plot(ax1, X, Y, 'b');
    hold on

    if(f<0.1)
        gname = ['disp_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['disp_volt_f', num2str(f*1000), '_O0', num2str(Omega*10)];
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
    close(fig5);
end
% ...........................................................

%% plot velo vs volt
if velo_volt == 1
    Nss = round(0.9*Ndt);

    X = Qvelo(Nss:Ndt);
    Y = Qvolt(Nss:Ndt);

    fig6 = figure('NumberTitle','off');
    ax1 = axes('Position',[0.2 0.2 0.7 0.7]);
    fh1 = plot(ax1, X, Y, 'b');
    hold on

    if(f<0.1)
        gname = ['velo_volt_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['velo_volt_f', num2str(f*1000), '_O0', num2str(Omega*10)];
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
    close(fig6);
end
% ...........................................................

%% plot atractor
if attractor == 1
    fig7 = figure('NumberTitle','off');

    X1 = Qdisp(Nss:Ndt);
    Y1 = Qvelo(Nss:Ndt);
    Z1 = Qvolt(Nss:Ndt);

    fh1 = plot3(X1, Y1, Z1, 'b');
    hold on

    if(f<0.1)
        gname = ['3d_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname = ['3d_f', num2str(f*1000), '_O0', num2str(Omega*10)];
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
    close(fig7);
end
% ...........................................................

%% plot histograms
if histogram == 1
    %Nbins = round(sqrt(Ndt));
    Nbins = 64;
    Nksd  = 32*Nbins;

    Ndt = length(Qvolt);
    Nss = round(0.9*Ndt);
    Pot = lambda.*Qvolt(Nss:Ndt).*Qvolt(Nss:Ndt);

    [Pot_bins,Pot_freq] = randvar_pdf(Pot,Nbins);
    [Pot_ksd ,Pot_supp] = randvar_ksd(Pot,Nksd);

    gtitle    = ' ';
    xlab      = ' power';
    ylab      = ' probability density';
    xmin      = 0.0;
    xmax      = max(Pot_bins);
    ymin      = 0.0;
    ymax      = 'auto';

    if(f<0.1)
        gname     = ['pdf_pot_f0', num2str(f*1000), '_O0', num2str(Omega*10)];
    else
        gname     = ['pdf_pot_f', num2str(f*1000), '_O0', num2str(Omega*10)];
    end

    linearity = 'nonlinear';
    flag      = 'png';
    fig8      = graph_bar_curve1(Pot_bins,Pot_freq,...
                              Pot_supp,Pot_ksd,gtitle,...
                              xlab,ylab,xmin,xmax,ymin,ymax,gname,linearity,flag);

    close(fig8);
end
end
%%
toc
% -----------------------------------------------------------