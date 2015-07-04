clear all; close all;
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
monkey='Quincy';
area='LIP'; %'PITd';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

load(cell_file);


all_cell_results=cell(length(cell_str),2);
%%
for cell_no=1:length(cell_str)
    if ~isempty(cell_str(cell_no).attention.mat);
        %[results_single]=make_GLM_fun(cell_no,monkey,area,'fixed_points');
        [results_single]=make_GLM_fun(cell_no,monkey,area,'time_course');
        close(gcf)
        all_cell_results(cell_no,:)=results_single;
    end
end


%%
