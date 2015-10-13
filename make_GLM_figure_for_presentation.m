clc; clear all; close all;
area='PITd';
monkey='Quincy';
cell_no=29;
use_BRTs=0;
save_figures=1;

figure_dir='/Freiwald/ppolosecki/harbor';
filename=[area '_' monkey '_cell_' num2str(cell_no)];

    true_contrasts=[1 1;
                     1 2;
                     1 10;
                     2 8;
                     2 5;
                     1 7;
                     2 11;
                     2 7];
                 
% true_contrasts=[1 1;
%                 1 2;
%                 1 10;
%                 2 8;
%                 2 5;
%                 1 7;
%                 2 7];

[results]=make_GLM_fun(cell_no,monkey,area,'time_course',1,1);
close all


f=plot_GLM_contrasts_for_presentation(results,true_contrasts,use_BRTs);


if save_figures
    for i=1:length(f)
        plot2svg(fullfile(figure_dir,[filename '_part_' num2str(i) '.svg']),f(i));
        %saveas(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.fig']));
        print(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.png']),'-dpng','-r300')
    end
    
end