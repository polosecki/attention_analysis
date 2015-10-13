close all; clear all
monkeys = {'Quincy','Michel'};
areas = {'PITd','LIP'};
save_fig = false;true;


addpath('..')
addpath('../export_fig/')
addpath('../plot2svg/')

for mm=1:2
    for aa=1:2
        load([monkeys{mm} '_' areas{aa} '_performance.mat'])
        t=[repmat({monkeys{mm}},size(behavior_performance))];
        [behavior_performance.monkey] = t{:};
        t=[repmat({areas{aa}},size(behavior_performance))];
        [behavior_performance.area] = t{:};
        
        if ~exist('total_behav','var')
            total_behav = behavior_performance;
        else
            total_behav = [total_behav behavior_performance];
        end
        clear behavior_performance
    end
end

%%
close all
nbins=103; % 29;
cmap_used='bone';'hot';
zlim=0.11%01;
do_smooth=true;
if do_smooth
    K_size=round(nbins/10);
    fwhm=min(5,K_size/2);
    [X,Y]=meshgrid(linspace(-10,10,K_size),linspace(-10,10,K_size));
    sigma=sqrt(log(2))*fwhm;
    k=exp(-(X.^2+Y.^2)/(sigma^2));
    k=k/sum(k(:));
    %N=conv2(N,k,'same'); clear k
end

f=figure;
subplot(2,2,1)
[RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',1);
%edges=linspace(min(RT_abs),max(RT_abs),nbins)';
edges=linspace(0,max(RT_abs),nbins)';

[N,bc]=hist3([RT_abs BRT_abs],'edges',{edges edges});
N=N/sum(N(:))*100; if do_smooth; N=conv2(N,k,'same');end
axis equal
ah=imagesc(N,[0 zlim]); 
colormap(cmap_used)
%colorbar
yt=get(gca,'YTick');
xt=get(gca,'YTick');
t=bc{2};
set(gca,'YTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(yt)),'UniformOutput',false));
t=bc{1};
set(gca,'XTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(xt)),'UniformOutput',false));
hold on
xlim([.5 size(N,1)+.5]);
ylim([.5 size(N,1)+.5]);
xl=xlim;
yl=ylim;
axis square
lh=line(xl,yl,'color','b','LineWidth',3);
xlabel('Cued Surface t_{PME\_OT} (s)')
ylabel('RT(s)')
title('M1')

subplot(2,2,3)
edges=linspace(min(nRT_abs),max(nRT_abs),nbins)';
[N,bc]=hist3([nRT_abs nBRT_abs],'edges',{edges edges});
N=N/sum(N(:))*100; if do_smooth; N=conv2(N,k,'same'); end 
axis equal
ah=imagesc(N,[0 zlim]); 
colormap(cmap_used)
%colorbar
yt=get(gca,'YTick');
xt=get(gca,'YTick');
t=bc{2};
set(gca,'YTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(yt)),'UniformOutput',false));
t=bc{1};
set(gca,'XTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(xt)),'UniformOutput',false));
hold on
xlim([.5 size(N,1)+.5]);
ylim([.5 size(N,1)+.5]);
xl=xlim;
yl=ylim;
axis square
lh=line(xl,yl,'color','r','LineWidth',3);
xlabel('Distractor Surface t_{PME\_OT} (s)')
ylabel('RT(s)')


subplot(2,2,2)
[RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',2);
edges=linspace(min(RT_abs),max(RT_abs),nbins)';
[N,bc]=hist3([RT_abs BRT_abs],'edges',{edges edges});
N=N/sum(N(:))*100; if do_smooth; N=conv2(N,k,'same'); end 
axis equal
ah=imagesc(N,[0 zlim]); 
colormap(cmap_used)
p=get(gca,'position'); % save position
c_h=colorbar;
set(gca,'position',p); % restore position
xlabel(c_h,'%')
%colorbar
yt=get(gca,'YTick');
xt=get(gca,'YTick');
t=bc{2};
set(gca,'YTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(yt)),'UniformOutput',false));
t=bc{1};
set(gca,'XTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(xt)),'UniformOutput',false));
hold on
xlim([.5 size(N,1)+.5]);
ylim([.5 size(N,1)+.5]);
xl=xlim;
yl=ylim;
axis square
lh=line(xl,yl,'color','b','LineWidth',3);
xlabel('Cued Surface t_{PME\_OT} (s)')
ylabel('RT(s)')
title('M2')

subplot(2,2,4)
edges=linspace(min(nRT_abs),max(nRT_abs),nbins)';
[N,bc]=hist3([nRT_abs nBRT_abs],'edges',{edges edges});
N=N/sum(N(:))*100; if do_smooth; N=conv2(N,k,'same'); end 
axis equal
ah=imagesc(N,[0 zlim]); 
colormap(cmap_used)
%colorbar
yt=get(gca,'YTick');
xt=get(gca,'YTick');
t=bc{2};
set(gca,'YTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(yt)),'UniformOutput',false));
t=bc{1};
set(gca,'XTickLabel',cellfun(@(x) sprintf('%1.1f',x),num2cell(t(xt)),'UniformOutput',false));
hold on
xlim([.5 size(N,1)+.5]);
ylim([.5 size(N,1)+.5]);
xl=xlim;
yl=ylim;
axis square
lh=line(xl,yl,'color','r','LineWidth',3);
xlabel('Distractor Surface t_{PME\_OT} (s)')
ylabel('RT(s)')

if save_fig
    export_fig /Freiwald/ppolosecki/harbor/join_behav_distributions_figure.eps -eps -transparent
end