function save_many_GLMs(monkey,area,is_highres_GLM)
if nargin<3
    is_highres_GLM=0;
end
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
%monkey='Quincy';
%area='PITd'; %'PITd';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
%results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

load(cell_file);
addpath('..')

%all_cell_results=cell(length(cell_str),2);
%%
all_cell_results={};
for cell_no=1:length(cell_str)
    if ~ismember(cell_no,bad_files(monkey,area))
        %[results_single]=make_GLM_fun(cell_no,monkey,area,'fixed_points');
        %[results_single]=make_GLM_fun(cell_no,monkey,area,'time_course',1,1);
        if is_highres_GLM
        [results_single]=make_GLM_fun(cell_no,monkey,area,'betas_for_pca',1,1);
        else
        [results_single]=make_GLM_fun(cell_no,monkey,area,'time_course',1,1); 
        end        
        close(gcf)
        drawnow
        for k=1:length(results_single)
            results_single{k}.area=area;
            results_single{k}.monkey=monkey;
            results_single{k}.cell_no=cell_no;
        end
        all_cell_results=[all_cell_results;results_single];        
    end
end

%save([area '_' monkey '_GLMs.mat'],'all_cell_results')
if is_highres_GLM
save([area '_' monkey '_highres_GLMs.mat'],'all_cell_results','-v7.3')
else
save([area '_' monkey '_GLMs.mat'],'all_cell_results')
end
%%
