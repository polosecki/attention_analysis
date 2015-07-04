close all
area='LIP'; %'PITd';
monkey='Michel';
cell_no=28;%29;%07;
figure_dir='/home/ppolosecki/figures';
filename=[monkey '_' area '_' sprintf('%03.0f',cell_no) '_GLM'];
%[results]=make_GLM_fun(cell_no,monkey,area,'fixed_points');

[results]=make_GLM_fun(cell_no,monkey,area,'time_course');
%addpath(fullfile('.','plot2svg'))
%plot2svg(fullfile(figure_dir,[filename '.svg']),gcf);
%saveas(gcf,fullfile(figure_dir,[filename '.fig']));


 %contrasts_plotted={logical([1 1 0 0 0 0 1 0 0]);
 %logical([0 0 0 0 1 0 1 1 0])};
%  contrasts_plotted={logical([0 0 0 0 0 0 0 0 0]);
%  logical([0 0 0 0 0 0 0 0 0])};
true_contrasts=[1 1;
 1 2;
 2 8;
 2 5;
 2 7
 1 7;];
 
close all
f=plot_GLM_contrasts_for_presentation(results,true_contrasts);
for i=1:length(f)
plot2svg(fullfile(figure_dir,[filename '_part_' num2str(i) '.svg']),f(i));
saveas(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.fig']));
print(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.png']),'-dpng','-r300')
end

%plot2svg(['Michel_LIP_surface_cells' '.svg'],gcf);
%%

if false
    [results]=make_GLM_fun(cell_no,monkey,area,'fixed_points');
    
    contrasts_plotted={logical([1 1 0 0 0 0 1 0 0]);
        logical([0 0 0 0 1 0 1 1 0])};
    colors_used=distinguishable_colors(sum(contrasts_plotted{1}+sum(contrasts_plotted{2})));
    f=figure;
    for mat_used=1:2
        subplot(1,2,mat_used)
        bar_mat=[results{mat_used}.GLM1.ces(contrasts_plotted{1},:); results{mat_used}.GLM2.ces(contrasts_plotted{2},:)];
        ax=bar(bar_mat');
        for i=1:length(ax)
            set(ax(i),'FaceColor',colors_used(i,:))
        end
        set(gca,'XTickLabel',results{mat_used}.time)
        if mat_used==1
            legend([results{1}.GLM1.contrast.name(contrasts_plotted{1}); results{1}.GLM2.contrast.name(contrasts_plotted{2})]);
        end
    end
end