close all; clear all

group_histograms_attention

%load Quincy_LIP
%good_files=logical(ones(size(all_cell_results,1),1));
%good_files([2 12:14 17 18 24:26])=0;

%%temp=all_cell_results;

%load Michel_LIP
%good_files2=zeros(size(all_cell_results,1),1);
%good_files2([6 11 13 16 17 18 20 22:30])=1; %all visual files

all_cell_results=[temp;all_cell_results];
good_files=[good_files;good_files2];
clear good_files2 temp norm_group
scratch_population_analysis_pp

scratch_population_analysis_pp;