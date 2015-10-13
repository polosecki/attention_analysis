clear all; close all
%areas={'PITd','LIP'}
%monkeys={'Quincy','Michel'}
%for mm=1:length(monkey)
%    for aa=1:length(area);
area='PITd';%areas{aa};
monkey='Quincy';%'Michel';%'Quincy';%'both';%'Michel';%monkeys{mm};%'Quincy';
noise_model='poisson';
remake_GLMs=0;
save_fig=1;
use_highres_data=1; % TIME-CONSUMING!!!
use_BRTs=1;
use_BRTs_in_plot=0;

addpath('..');

if use_highres_data
    file_suffix= '_highres_GLMs.mat';
else
    file_suffix='_GLMs.mat';
end

if ~strcmp(monkey,'both')
    if ~exist([area '_' monkey '_' noise_model file_suffix],'file') | remake_GLMs
        save_many_GLMs(monkey,area,use_highres_data,noise_model)
    end
    load([area '_' monkey '_' noise_model file_suffix])
else
    load([area '_Quincy_' noise_model file_suffix]);
    temp=all_cell_results;
    load([area '_Michel_' noise_model file_suffix]);
    all_cell_results=[temp; all_cell_results];
    clear temp
end


if use_BRTs
    end_indx=size(all_cell_results,2);
else
    end_indx=2;
end

temp=[];
for k=1:end_indx
    temp=[temp [all_cell_results{:,k}]'];
end
%%
%Calculate normalization factor:
%OBSOLETE: refer to ../scratch_population_analysis_pp.m for original
%version
%
% aa=temp(:,1);
% temp2=[aa.GLM]; clear aa
% temp2=temp2(1:length(temp(1).GLM):end); %isolate GLM1
% %temp2=[temp(:,1).GLM(1)]; % GLM used
% norm_factor=nan(length(temp2),1);
% baseline_indexes=(temp(1,1).time >=-1.4) & (temp(1,1).time <=-1.1);
% for i=1:length(temp2)
% norm_factor(i)=mean(temp2(i).ces(9,baseline_indexes));
% end


for mat_used=1:end_indx
    norm_group{mat_used}.time=calculate_time_axis_group_GLM({temp(:,mat_used).time}',mat_used);
    
    for k=1:length(temp(1,1).GLM)
        norm_group{mat_used}.GLM(k).ces=nan(size(temp(1,1).GLM(k).ces,1),length(norm_group{mat_used}.time));
        norm_group{mat_used}.GLM(k).ces_std=nan(size(temp(1,1).GLM(k).ces,1),length(norm_group{mat_used}.time));
        
        
        for c=1:size(temp(1,1).GLM(k).ces,1)
            allces=nan(size(temp,1),length(norm_group{mat_used}.time));
            for i=1:size(temp,1)
                if mat_used==1
                    allces(i,1:length(temp(i,mat_used).time))=temp(i,mat_used).GLM(k).ces(c,:);
                elseif mat_used==2
                    allces(i,end-(length(temp(i,mat_used).time)-1):end)=temp(i,mat_used).GLM(k).ces(c,:);
                elseif mat_used==3
                    [~,bin_start]=min(abs(norm_group{mat_used}.time-temp(i,mat_used).time(1)));
                    bin_end=bin_start+length(temp(i,mat_used).time)-1;
                    if bin_end>(size(allces,2)+1)
                        error('Something fishy going on')
                    elseif bin_end==(size(allces,2)+1)
                     shortened_ces=temp(i,mat_used).GLM(k).ces(c,:);
                     shortened_ces=shortened_ces(1:end-1);
                     allces(i,bin_start:bin_end-1)=shortened_ces;
                    else %normal case
                    allces(i,bin_start:bin_end)=temp(i,mat_used).GLM(k).ces(c,:);
                    end
                    
                end
            end
            norm_group{mat_used}.GLM(k).ces(c,:)=nanmean(allces,1);
            norm_group{mat_used}.GLM(k).ces_std(c,:)=nanstd(allces,1)./sqrt(sum(~isnan(allces),1));
        end
        norm_group{mat_used}.GLM(k).contrast.name=all_cell_results{1,mat_used}.GLM(k).contrast.name;
    end
   
end

make_common_plot=1;
for_paper=1;
if make_common_plot
    contrasts_plotted={logical([0 0 0 0 0 0 0 0 0 0]);
        logical([0 0 0 0 0 0 0 0 0]);
        logical([1 1 1 0 1 0 0 1 1 0 1 0 1 0 0 0])};
%     contrasts_plotted={logical([1 0 0 0]);
%         logical([0 1 0 0]);
%         logical([1 1 0 0 0 0 0 1 1 0 1 0 1 0 0 0])};
    if ~for_paper
        h=plot_GLM_contrasts(norm_group,contrasts_plotted,use_BRTs_in_plot);
    else
%        contrasts_plotted={false(1,length(norm_group{2}.GLM(1).contrast.name)),
%    false(1,length(norm_group{2}.GLM(2).contrast.name))
%    false(1,length(norm_group{2}.GLM(3).contrast.name))};
%contrasts_plotted{3}([3,5,13])=true;
        h=plot_GLM_contrasts_for_paper(norm_group,contrasts_plotted,use_BRTs_in_plot);
    end
    [axh,labelh]=suplabel([area '_' monkey],'t');
    set(labelh,'Interpreter','none');
    drawnow
    if save_fig
        figure_dir='/Freiwald/ppolosecki/harbor/';
        filename=[area '_' monkey '_group_GLM'];
        addpath(fullfile('..','plot2svg'));
        plot2svg(fullfile(figure_dir,[filename '.svg']),gcf);
    end
else
    figure_dir='/Freiwald/ppolosecki/harbor/';
    filename=[area '_' monkey '_group_GLM'];
    true_contrasts=[1 1;
                     1 2;
                     1 10;
                     2 8;
                     2 5;
                     1 7;
                     2 11;
                     2 7];                                                 
%     true_contrasts=[1 1;
%         1 2;
%         2 8;
%         2 5;
%         2 7
%         1 7;];
    close all
    f=plot_GLM_contrasts_for_presentation(norm_group,true_contrasts,use_BRTs_in_plot);
    for i=1:length(f)
        plot2svg(fullfile(figure_dir,[filename '_part_' num2str(i) '.svg']),f(i));
        %saveas(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.fig']));
        print(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.png']),'-dpng','-r300')
    end
end

%    end
%end



