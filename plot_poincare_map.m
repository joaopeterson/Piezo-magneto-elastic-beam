
% -----------------------------------------------------------------
%  plot_poincare_map.m
%
%  This functions plots the Poincare map of dynamics system.
%
%  input:
%  x1     - x data vector 1
%  y1     - y data vector 1
%  xmap   - x Poincare map
%  ymap   - y Poincare map
%  gtitle - graph title
%  leg1   - legend 1
%  leg2   - legend 2
%  leg3   - legend 3
%  xlab   - x axis label
%  ylab   - y axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Fev 10, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = plot_poincare_map(x1,y1,xmap,ymap,gtitle,...
                                 xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)
    
    % check number of arguments
    if nargin < 12
        error('Too few inputs.')
    elseif nargin > 13
        error('Too many inputs.')
    elseif nargin == 12
        flag = 'none';
    end

    % check arguments
    if length(x1) ~= length(y1)
        error('x1 and y1 vectors must be same length')
    end
    
    if length(xmap) ~= length(ymap)
        error('xmap and ymap vectors must be same length')
    end
    
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = plot(x1,y1,'-b');
    hold all
    fh2 = plot(xmap,ymap,'.');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    axis equal
    
    if ( strcmp(xmin,'auto') || strcmp(xmax,'auto') )
        xlim('auto');
    else
        xlim([xmin xmax]);
    end
    
    if ( strcmp(ymin,'auto') || strcmp(ymax,'auto') )
        ylim('auto');
    else
        ylim([ymin ymax]);
    end
    
    set(fh1,'LineWidth',0.3);
    set(fh2,'MarkerSize',4);
    set(fh2,'MarkerEdgeColor','r');
    labX = xlabel(xlab,'FontSize',22,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',22,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
    hold off
    
	title(gtitle,'FontSize',22,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        graph_fixPSlinestyle(gname,gname);
        
    elseif (strcmp(flag,'png'))
        saveas(gcf, gname, 'png');
    end

return
% -----------------------------------------------------------------
