function [cell_ces, cell_sig, param_names]=calculate_stats_from_cell_list(monkey,area)

cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
load(cell_file)
good_files=true(length(cell_str),1);
% if strcmp(monkey,'Quincy') && strcmp(area,'PITd')
%     good_files([4 8 22 23 24 32 33 42:46 52])=false;
% elseif strcmp(monkey,'Quincy') && strcmp(area,'LIP')
%     good_files([2 12:14 17 26])=false;
% elseif strcmp(monkey,'Michel') && strcmp(area,'LIP')
%     good_files([1 2 4 5 21 31 32 33])=0; % cell 1 removed beause it is a copy of cell 12
% elseif strcmp(monkey,'Michel') && strcmp(area,'PITd')
%     good_files([6 10:12 16 19 32 38 40 48])=0;
%end
good_files(bad_files(monkey,area))=false;

for cell_no=1:length(cell_str)
    if good_files(cell_no)
        [results]=make_GLM_fun(cell_no,monkey,area,'fixed_points',1,0);
        %[results]=make_GLM_fun(cell_no,monkey,area,'betas_for_pca',1,0);
        if ~isempty(cell_str(cell_no).MGS_file.mat)
            wd=pwd;
            cd ('../../mgs_analysis')
            %[results_mgs]=make_GLM_mgs_fun(cell_no,monkey,area,'betas_for_pca',1);
            [results_mgs]=make_GLM_mgs_fun(cell_no,monkey,area,'fixed_points',1,1);
            cd(wd)
            %cell_results_collection_mgs(cell_no,1:2)=[results_mgs{1}.GLM results_mgs{2}.GLM];
            cell_results_collection_mgs(cell_no,1:3)=[results_mgs{1}.GLM results_mgs{2}(1).GLM results_mgs{2}(2).GLM];
        end
        cell_results_collection(cell_no,1:3)=results{2}.GLM;
        firing_mean_std(cell_no,:)=[results{2}.mean_activity results{2}.std_activity];
        disp(['cell done: ' num2str(cell_no)])

    end
end

%%
%Surf: cell_results_collection(1,1).contrast.name(1)
%Target:cell_results_collection(1,1).contrast.name(2)
%Attn(Targets out of RF): cell_results_collection(1,1).contrast.name(4)
%Attn/Target inter: cell_results_collection(1,1).contrast.name(7)
%Saccade: cell_results_collection(1,2).contrast.name(5)
%Saccade/attention inter: cell_results_collection(1,2).contrast.name(7)
%Guide to extract appropiate effect sizes:
%cell_results_collection(cell_no,GLM_no).ces(contrast_number,time_point)


%parameters of interest
% params=[1 4 1; %attn: GLM1,4thcontrast,time_point1
%         1 2 1; %target
%         1 7 1; %attn/target inter
%         2 5 2; %saccade
%         2 7 2; %saccade/attn inter
%         1 1 1; %surface
%         1 10 1; %surface/target inter
%         2 11 1]; %saccade/surface interaction
    
params=[3 3 1; %attn: GLM3,3rdcontrast,time_point1
        3 2 1; %target
        3 11 1; %attn/target inter
        3 5 2; %saccade
        3 13 2; %saccade/attn inter
        3 1 1; %surface
        3 8 1; %surface/target inter
        3 9 1]; %saccade/surface interaction     
param_names={'attention';
             'target';
             'attn/target inter';
             'saccade';
             'saccade/attn inter';
             'surface';
             'surf/targ inter';
             'saccade/surf inter';
             'mean firing';
             'std firing';
             'mgs visual epoch';
             'mgs memory poch';
             'mgs saccade epoch';
             'mgs mem-sacc diff'};

cell_ces=nan(size(cell_results_collection,1),size(params,1)+5);
cell_sig=nan(size(cell_results_collection,1),size(params,1)+5);
for cell_no=1:size(cell_results_collection,1)
    if ~isempty(cell_results_collection(cell_no,1).ces)
        for par_num=1:size(params,1)
            cell_ces(cell_no,par_num)=cell_results_collection(cell_no,params(par_num,1)).ces(params(par_num,2),params(par_num,3));
            cell_sig(cell_no,par_num)=cell_results_collection(cell_no,params(par_num,1)).Fsig(params(par_num,2),params(par_num,3));
        end
            cell_ces(cell_no,par_num+[1:2])=[firing_mean_std(cell_no,1) firing_mean_std(cell_no,2)/firing_mean_std(cell_no,1)];
            cell_sig(cell_no,par_num+[1:2])=[0.01 0.01];

    end
    if ~isempty(cell_results_collection_mgs(cell_no,1).ces)
        cell_ces(cell_no,par_num+3)=cell_results_collection_mgs(cell_no,1).ces;
        cell_sig(cell_no,par_num+3)=cell_results_collection_mgs(cell_no,1).Fsig;
        cell_ces(cell_no,par_num+[4:5])=cell_results_collection_mgs(cell_no,2).ces;
        cell_sig(cell_no,par_num+[4:5])=cell_results_collection_mgs(cell_no,2).Fsig;
        cell_ces(cell_no,par_num+6)=cell_results_collection_mgs(cell_no,3).ces;
        cell_sig(cell_no,par_num+6)=cell_results_collection_mgs(cell_no,3).Fsig;
    end
end
%sendmail('9178422510@txt.att.net', 'Function analyses report', 'Loop completed');