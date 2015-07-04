function fh = plot_2d_distribution_with_marginals(N,bc,options)

if isfield(options,'smooth_N') && options.smooth_N==true
    K_size=round(size(N,1)/10);
    fwhm=min(5,K_size/2);
    [X,Y]=meshgrid(linspace(-10,10,K_size),linspace(-10,10,K_size));
    sigma=sqrt(log(2))*fwhm;
    k=exp(-(X.^2+Y.^2)/(sigma^2));
    k=k/sum(k(:));
    N=conv2(N,k,'same'); clear k
end

fh = figure;
set(fh, 'menubar', 'none');
%h(1)=subplot(3,3,[1 4]);
h(1)=subplot(3,3,[4 7]);

%bar(bc{1},sum(N,2))
%set(get(h(1),'children'),'FaceColor',[.5 .5 .5])
plot(bc{1},sum(N,2),'k-')
xlim([bc{1}(1) bc{1}(end)])
set(h(1),'view',[-90 -90]) %azimuth, elevation
if isfield(options,'ylim')
    ylim(options.ylim)
end
%h(2)=subplot(3,3,[8 9]);
h(2)=subplot(3,3,[2 3]);

%bar(bc{1},sum(N,1))
%set(get(h(2),'children'),'FaceColor',[.5 .5 .5])
plot(bc{1},sum(N,1),'k-')
xlim([bc{1}(1) bc{1}(end)])

%set(h(2),'view',[0 -90]) %azimuth, elevation
set(h(2),'view',[0 90]) %azimuth, elevation

if isfield(options,'xlim')
    ylim(options.xlim)
end

%h(3)=subplot(3,3,[2 3 5 6]);
h(3)=subplot(3,3,[5 6 8 9]);

if isfield(options,'zlim')
    ah=imagesc(N,options.zlim);
else
    ah=imagesc(N);
end
if isfield(options,'title')
    title(h(3),options.title)
end
if isfield(options,'cmap')
    colormap(h(3),options.cmap);
end
axis square
[~,b]=min(abs(bc{2}));
xl=xlim;
yl=ylim;
lh=line([b b],yl,'color','w','LineWidth',1);
lh=line(xl,[b b],'color','w','LineWidth',1);
lh=line(xl,yl,'color','w','LineWidth',1);
p=get(h(3),'position'); % save position
c_h=colorbar;
xlabel(c_h,'%')
set(h(3),'position',p); % restore position
drawnow
axis off

if isfield(options,'diag_hist')
%    h(4)=subplot(3,3,7);
    h(4)=subplot(3,3,1);
    [nn,xc]=hist(options.diag_hist,size(N,1));
    if isfield(options,'smooth_N') && options.smooth_N==true
    K_size=round(size(N,1)/10);
    X=linspace(-10,10,K_size);
    fwhm=max(5,K_size/2);
    sigma=sqrt(log(2))*fwhm;
    k=exp(-(X.^2)/(sigma^2));
    k=k/sum(k(:));
    nn=conv(nn,k,'same');
    end
    plot(xc,nn/sum(nn)*100,'k')
    %set(get(h(4),'children'),'FaceColor',[.5 .5 .5])
    line([0 0], ylim, 'color','k')
end

tp=get(h(3),'Position');
for i=1:2
    cp=get(h(i),'Position');
    %tda=get(h(3),'PlotBoxAspectRatio')
    %cda=get(h(i),'PlotBoxAspectRatio')
    
    if i==1
        cp(4)=tp(4);
        cp(1)=cp(1)+0.109;
        a=get(h(i),'ytick');
        set(h(i),'ytick',a(2:end));
    elseif i==2
        cp(1)=cp(1)+.055;
        cp(2)=tp(2)+tp(4)+0.01;
        cp(3)=cp(3)-.11;
        cp(4)=cp(4)+.0562;
        
        %     cp(3)=tp(3);
    end
    set(h(i),'Position',cp);
    set(h(i),'xlim',[bc{1}(1) bc{1}(end)]);
    set(get(h(i),'xlabel'),'string','time(s)')
    if i==1
        set(h(i),'YAxisLocation','right')
        if isfield(options,'ylabel')
            set(get(h(i),'xlabel'),'string',options.ylabel)
        else
            set(get(h(i),'xlabel'),'string','time(s)')
        end
            set(get(h(i),'xlabel'),'rotation',90)
    end
    
    if i==2
        set(h(i),'YAxisLocation','right')
        set(h(i),'XAxisLocation','top')

        if isfield(options,'xlabel')
            set(get(h(i),'xlabel'),'string',options.xlabel)
        else
            set(get(h(i),'xlabel'),'string','time(s)')
        end
    end
    set(get(h(i),'ylabel'),'string','%');
    
end

if length(h)>=4
tp=get(h(4),'Position');
cp=get(h(1),'Position');
tp(1)=cp(1);
cp=get(h(2),'Position');
tp(4)=cp(4);
if isfield(options,'rotate_diagonal') && options.rotate_diagonal==true
tp(2)=tp(2)-0.05;
tp(1)=tp(1)+0.05;
tp(2)=tp(2)-0.10;
end
set(h(4),'Position',tp);
%set(h(4),'view',[-80 90])
end

if isfield(options,'no_axis') && options.no_axis==true
for i=1:length(h)
    axis(h(i),'off')
end
end