clc;clear all; close all
monkey='Michel';
area='PITd';
cell_no=35;
[RT_out conditions_plotted_names conditions_plotted_for_computation]=measure_RTs_attention_for_cell(cell_no,monkey,area);
%%
y=[];phi=[];brt=[];
for i=1:length(RT_out)
    y=[y;RT_out{i}];
    phi=[phi;conditions_plotted_names(i,1)*ones(size(RT_out{i}))];
    brt=[brt; conditions_plotted_names(i,2)*ones(size(RT_out{i}))];
end

[p,table,stats]=anovan(y,{phi brt},'model','full');
 
%c=multcompare(stats)

%% 
figure;
load('CoolWarmMap')
scale_semiwidth=2*std(cellfun(@median,RT_out));
scale_center=mean(cellfun(@median,RT_out));
imagesc(reshape(cellfun(@median,RT_out),[4 4])',[scale_center-scale_semiwidth scale_center+scale_semiwidth])
colormap(cmap)
set(gca,'xtick',1:4)
set(gca,'xticklabel',unique(conditions_plotted_names(:,1)))
set(gca,'ytick',1:4)
set(gca,'yticklabel',unique(conditions_plotted_names(:,1)))
ylabel('Attended Surface')
xlabel('BRT direction')
title('Reaction times')

