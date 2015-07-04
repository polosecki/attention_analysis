clear all; close all
monkey='Quincy';
area='PITd';
cell_no=29;

overwrite_grand_psth=1;
figure_dir=[];


%
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
base_dir=fullfile('/Freiwald/ppolosecki/lspace/',lower(monkey));
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);
load(results_file)
RF_pos=results_data(cell_no).RF_pos;


[closest_surf_phi] = attention_analysis_func_pp(cell_file,cell_no,base_dir,overwrite_grand_psth,figure_dir,RF_pos);

