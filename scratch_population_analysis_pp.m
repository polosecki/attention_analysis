% %Quincy PITd good cells
% good_files=logical(ones(size(all_cell_results,1),1));
% good_files([4 8 22 23 24 32 33 42:46 52])=0;


 %Quincy LIP good cells
good_files=logical(ones(size(all_cell_results,1),1));
good_files([2 12:14 17 26])=0;


% %Michel LIP good cells
%good_files=logical(ones(size(all_cell_results,1),1));
%good_files([2 4 5 21 31 32 33])=0;

% good_files([6 11 13 16 17 18 20 22:30])=1; %all visual files
%good_files([6 13 16 17 22 23 30])=1; %all surface-like files
%good_files([11 18 24:29])=1; %all target-like files

% %Michel PITd
%good_files=logical(ones(size(all_cell_results,1),1));
%good_files([6 10:12 16 19 32 38 40 48])=0;


temp=[[all_cell_results{good_files,1}]' [all_cell_results{good_files,2}]'];

%Calculate normalization factor:
aa=temp(:,1);
temp2=[aa.GLM]; clear aa
temp2=temp2(1:length(temp(1).GLM):end); %isolate GLM1
%temp2=[temp(:,1).GLM(1)]; % GLM used
norm_factor=nan(length(temp2),1);
baseline_indexes=(temp(1,1).time >=-1.4) & (temp(1,1).time <=-1.1);
for i=1:length(temp2)
norm_factor(i)=mean(temp2(i).ces(9,baseline_indexes)); 
end

norm_group=all_cell_results(1,:);%results_single;
for mat_used=1:2
    all_times={temp(:,mat_used).time}';
    [vm,Im]=max(cellfun(@length,all_times));
    norm_group{mat_used}.time=all_times{Im};
    
    norm_group{mat_used}.GLM(1).ces=nan(size(norm_group{mat_used}.GLM(1).ces,1),length(norm_group{mat_used}.time));
    norm_group{mat_used}.GLM(1).ces_std=nan(size(norm_group{mat_used}.GLM(1).ces,1),length(norm_group{mat_used}.time));
    norm_group{mat_used}.GLM(2).ces=nan(size(norm_group{mat_used}.GLM(2).ces,1),length(norm_group{mat_used}.time));
    norm_group{mat_used}.GLM(2).ces_std=nan(size(norm_group{mat_used}.GLM(2).ces,1),length(norm_group{mat_used}.time));
    
    for c=1:size(norm_group{mat_used}.GLM(1).ces,1)
        allces=nan(size(temp,1),length(norm_group{mat_used}.time));
        for i=1:size(temp,1)
            if mat_used==1
                allces(i,1:length(temp(i,mat_used).time))=temp(i,mat_used).GLM(1).ces(c,:)./norm_factor(i);
            elseif mat_used==2
                allces(i,end-(length(temp(i,mat_used).time)-1):end)=temp(i,mat_used).GLM(1).ces(c,:)./norm_factor(i);
            end
        end
        norm_group{mat_used}.GLM(1).ces(c,:)=nanmean(allces,1);
        norm_group{mat_used}.GLM(1).ces_std(c,:)=nanstd(allces,1)./sqrt(sum(~isnan(allces),1));
    end
    for c=1:size(norm_group{mat_used}.GLM(2).ces,1)
        allces=nan(size(temp,1),length(norm_group{mat_used}.time));
        for i=1:size(temp,1)
            if mat_used==1
                allces(i,1:length(temp(i,mat_used).time))=temp(i,mat_used).GLM(2).ces(c,:)./norm_factor(i);
            elseif mat_used==2
                allces(i,end-(length(temp(i,mat_used).time)-1):end)=temp(i,mat_used).GLM(2).ces(c,:)./norm_factor(i);
            end
        end
        norm_group{mat_used}.GLM(2).ces(c,:)=nanmean(allces,1);
        norm_group{mat_used}.GLM(2).ces_std(c,:)=nanstd(allces,1)./sqrt(sum(~isnan(allces),1));
        
    end
end

make_common_plot=1;
if make_common_plot
    contrasts_plotted={logical([1 1 0 0 0 0 1 0 0 1]);
        logical([0 0 0 0 1 0 1 1 0])};
    plot_GLM_contrasts(norm_group,contrasts_plotted);
else
    figure_dir='/home/ppolosecki/figures';
    filename='group_GLM';
    true_contrasts=[1 1;
        1 2;
        2 8;
        2 5;
        2 7
        1 7;];
    close all
    f=plot_GLM_contrasts_for_presentation(norm_group,true_contrasts);
    for i=1:length(f)
        plot2svg(fullfile(figure_dir,[filename '_part_' num2str(i) '.svg']),f(i));
        saveas(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.fig']));
        print(f(i),fullfile(figure_dir,[filename '_part_' num2str(i) '.png']),'-dpng','-r300')
    end
end




