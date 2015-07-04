function [ax,hg1,hg2] = mzplotyy(x,y1,y2,col1,col2,ax);

create_axes = 0;
if ~exist('ax','var') | isempty(ax)
    create_axes = 1;
end

if create_axes
    cax = newplot;
    fig = get(cax,'Parent');
    ax(1) = cax;
    set(fig,'NextPlot','add')
else
    axes(ax(1));
    cax = ax(1);
    fig = get(cax,'Parent');
end

% Plot first plot
hg1 = plot(x,y1,'Color',col1);
set(ax(1),'Color','none');
set(ax(1),'Box','off');

if create_axes
    ax1hv = get(ax(1),'HandleVisibility');
    ax(2) = axes('HandleVisibility',ax1hv,'Units',get(ax(1),'Units'), ...
        'Position',get(ax(1),'Position'),'Parent',fig);
else
    axes(ax(2));
end
hg2 = plot(x,y2,'Color',col2);
set(ax(2),'XAxisLocation','top','YAxisLocation','right','Color','none',...
        'XGrid','off','YGrid','off','Box','off');

axes(ax(1));

% create DeleteProxy objects (an invisible text object in
% the first axes) so that the other axes will be deleted
% properly.
if create_axes
    DeleteProxy(1) = text('Parent',ax(1),'Visible','off',...
                          'Tag','MzPlotyyDeleteProxy',...
                          'HandleVisibility','off',...
            'DeleteFcn','try;delete(get(gcbo,''userdata''));end');
    DeleteProxy(2) = text('Parent',ax(2),'Visible','off',...
                           'Tag','MzPlotyyDeleteProxy',...
                           'HandleVisibility','off',...
            'DeleteFcn','try;delete(get(gcbo,''userdata''));end');
    set(DeleteProxy(1),'UserData',ax(2));
    set(DeleteProxy(2),'UserData',DeleteProxy(1));
end
