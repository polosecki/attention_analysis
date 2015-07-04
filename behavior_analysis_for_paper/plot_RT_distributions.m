% Load data

close all; clear all
monkeys = {'Quincy','Michel'};
areas = {'PITd','LIP'};
save_fig = false;


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
%% Make basic histograms of BRTs, RTs
n_bins=99;
bin_edges=linspace(0,6,n_bins+1);
bin_centers=(bin_edges(1:end-1)+bin_edges(2:end))/2;

[RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',3);
f=figure;
subplot(2,2,1)
[n,x] = hist(BRT_abs,bin_centers);
h=bar(x,n/sum(n)*100);
set(get(h,'children'),'FaceColor',[.5 .5 .5])
xlim([0 6])
title('t_{PME\_OT} distribution (cued surface)')
xlabel('Time from surface onset (s)')
ylabel('%')
ylim([0 5])
box off 

subplot(2,2,2)
[n,x] = hist(nBRT_abs,bin_centers);
h=bar(x,n/sum(n)*100);
set(get(h,'children'),'FaceColor',[.5 .5 .5])
xlim([0 6])
title('t_{PME\_OT} distribution (distractor surface)')
xlabel('Time from surface onset (s)')
ylabel('%')
ylim([0 5])
box off

subplot(2,2,3)
[RT_abs, ~, ~, ~] = get_distribution_from_total_behav(total_behav,'none',1);
[n,x] = hist(RT_abs,bin_centers);
h=bar(x,n/sum(n)*100);
set(get(h,'children'),'FaceColor',[.5 .5 .5])
xlim([0 6])
title('RT  distribution (M1)')
xlabel('Time from surface onset (s)')
ylabel('%')
ylim([0 5])
box off

subplot(2,2,4)
[RT_abs, ~, ~, ~] = get_distribution_from_total_behav(total_behav,'none',2);
[n,x] = hist(RT_abs,bin_centers);
h=bar(x,n/sum(n)*100);
set(get(h,'children'),'FaceColor',[.5 .5 .5])
xlim([0 6])
title('RT  distribution (M2)')
xlabel('Time from surface onset (s)')
ylabel('%')
ylim([0 5])
box off

if save_fig
    export_fig /Freiwald/ppolosecki/harbor/Raw_behav_distributions_figure.eps -eps -transparent
    pause(4)
    plot2svg('/Freiwald/ppolosecki/harbor/Raw_behav_distributions_figure.svg',f)
end

%% Make histograms of RTs condidioned on different BRT times
n_bins=40;
bin_edges=linspace(0,6,n_bins+1);
bin_centers=(bin_edges(1:end-1)+bin_edges(2:end))/2;
colors_used=distinguishable_colors(2);
f=figure;
set(f,'Position',get(0,'ScreenSize'));
hs=subplot(2,2,1);
box off
[RT_abs, ~, nRT_abs, ~] = get_distribution_from_total_behav(total_behav,'low',1);
[~, BRT_abs, ~, ~] = get_distribution_from_total_behav(total_behav,'none',1,false);

[n,x] = hist(RT_abs,bin_centers);
h1=bar(x,n/sum(n)*100);
set(get(h1,'children'),'FaceColor',colors_used(1,:),'EdgeColor',colors_used(1,:),'facealpha',0.75)
hold on
[n,x] = hist(nRT_abs,bin_centers);
h2=bar(x,n/sum(n)*100);
set(get(h2,'children'),'FaceColor',colors_used(2,:),'EdgeColor',colors_used(2,:),'facealpha',0.5)
xlim([0 6])
ylim([0 20])
xlabel('Time from surface onset (s)')
ylabel('%')
title('M1 (t_{PME\_OT} in first decile)')
yl=ylim;
vtx=[min(BRT_abs) min(BRT_abs) prctile(BRT_abs,10) prctile(BRT_abs,10)];
vty=[yl(1) yl(2) yl(2) yl(1)];
hp=patch(vtx,vty,[.5 .5 .5]);
set(hp,'facealpha',.5)
%line([prctile(BRT_abs,10) prctile(BRT_abs,10) ], ylim, 'Color', 'k');
%line([min(BRT_abs) min(BRT_abs)], ylim, 'Color', 'k');
legend_vals={'Using cued surface t_{PME\_OT}','Using cued surface t_{PME\_OT}','t_{PME\_OT} range'};
line_handles=[get(h1,'children'),get(h2,'children'),hp];
lgh = legend(line_handles,legend_vals);
set(lgh,'Box','off');

hs=subplot(2,2,2);
box off
[RT_abs, ~, nRT_abs, ~] = get_distribution_from_total_behav(total_behav,'low',2);
[~, BRT_abs, ~, ~] = get_distribution_from_total_behav(total_behav,'none',2,false);

[n,x] = hist(RT_abs,bin_centers);
h1=bar(x,n/sum(n)*100);
set(get(h1,'children'),'FaceColor',colors_used(1,:),'EdgeColor',colors_used(1,:),'facealpha',0.75)
hold on
[n,x] = hist(nRT_abs,bin_centers);
h2=bar(x,n/sum(n)*100);
set(get(h2,'children'),'FaceColor',colors_used(2,:),'EdgeColor',colors_used(2,:),'facealpha',0.5)
xlim([0 6])
ylim([0 20])
xlabel('Time from surface onset (s)')
ylabel('%')
title('M2 (t_{PME\_OT} in first decile)')
yl=ylim;
vtx=[min(BRT_abs) min(BRT_abs) prctile(BRT_abs,10) prctile(BRT_abs,10)];
vty=[yl(1) yl(2) yl(2) yl(1)];
hp=patch(vtx,vty,[.5 .5 .5]);
set(hp,'facealpha',.5)
%line([prctile(BRT_abs,10) prctile(BRT_abs,10) ], ylim, 'Color', 'k');
%line([min(BRT_abs) min(BRT_abs)], ylim, 'Color', 'k');
legend_vals={'Using cued surface t_{PME\_OT}','Using cued surface t_{PME\_OT}','t_{PME\_OT} range'};
line_handles=[get(h1,'children'),get(h2,'children'),hp];
lgh = legend(line_handles,legend_vals);
set(lgh,'Box','off')

hs=subplot(2,2,3);
box off
[RT_abs, ~, nRT_abs, ~] = get_distribution_from_total_behav(total_behav,'high',1);
[~, ~, ~, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',1,false);

[n,x] = hist(RT_abs,bin_centers);
h1=bar(x,n/sum(n)*100);
set(get(h1,'children'),'FaceColor',colors_used(1,:),'EdgeColor',colors_used(1,:),'facealpha',0.75)
hold on
[n,x] = hist(nRT_abs,bin_centers);
h2=bar(x,n/sum(n)*100);
set(get(h2,'children'),'FaceColor',colors_used(2,:),'EdgeColor',colors_used(2,:),'facealpha',0.5)
xlim([0 6])
ylim([0 20])
xlabel('Time from surface onset (s)')
ylabel('%')
title('M1 (t_{PME\_OT} in last decile)')
yl=ylim;
vtx=[prctile(BRT_abs,90) prctile(BRT_abs,90) max(BRT_abs) max(BRT_abs)];
vty=[yl(1) yl(2) yl(2) yl(1)];
hp=patch(vtx,vty,[.5 .5 .5]);
set(hp,'facealpha',.5)
%line([prctile(BRT_abs,90) prctile(BRT_abs,90) ], ylim, 'Color', 'k');
%line([max(BRT_abs) max(BRT_abs)], ylim, 'Color', 'k');
legend_vals={'Using cued surface t_{PME\_OT}','Using cued surface t_{PME\_OT}','t_{PME\_OT} range'};
line_handles=[get(h1,'children'),get(h2,'children'),hp];
lgh = legend(line_handles,legend_vals);
set(lgh,'Box','off','location','NorthWest')

hs=subplot(2,2,4);
[RT_abs, ~, nRT_abs, ~] = get_distribution_from_total_behav(total_behav,'high',2);
[~, ~, ~, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',2,false);
[n,x] = hist(RT_abs,bin_centers);
h1=bar(x,n/sum(n)*100);
set(get(h1,'children'),'FaceColor',colors_used(1,:),'EdgeColor',colors_used(1,:),'facealpha',0.75)
hold on
[n,x] = hist(nRT_abs,bin_centers);
h2=bar(x,n/sum(n)*100);
set(get(h2,'children'),'FaceColor',colors_used(2,:),'EdgeColor',colors_used(2,:),'facealpha',0.5)
xlim([0 6])
ylim([0 20])
xlabel('Time from surface onset (s)')
ylabel('%')
title('M2 (t_{PME\_OT} in last decile)')
%box off
yl=ylim;
vtx=[prctile(BRT_abs,90) prctile(BRT_abs,90) max(BRT_abs) max(BRT_abs)];
vty=[yl(1) yl(2) yl(2) yl(1)];
hp=patch(vtx,vty,[.5 .5 .5]);
set(hp,'facealpha',.5)
%line([prctile(BRT_abs,90) prctile(BRT_abs,90) ], ylim, 'Color', 'k');
%line([max(BRT_abs) max(BRT_abs)], ylim, 'Color', 'k');
legend_vals={'Using cued surface t_{PME\_OT}','Using cued surface t_{PME\_OT}','t_{PME\_OT} range'};
line_handles=[get(h1,'children'),get(h2,'children'),hp];
lgh = legend(line_handles,legend_vals);
set(lgh,'Box','off','location','NorthWest')

if save_fig
    export_fig /Freiwald/ppolosecki/harbor/RT_dist_figure.png -png -transparent
    pause(4)
    plot2svg('/Freiwald/ppolosecki/harbor/RT_dist_figure.svg',f)
end