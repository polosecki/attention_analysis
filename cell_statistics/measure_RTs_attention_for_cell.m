function [RT_out conditions_plotted_names conditions_plotted_for_computation]=measure_RTs_attention_for_cell(cell_no,monkey,area)
% Give single trial RTs for a given cell in a monkey and area. They are separated by condition
% INPUTS:
%cell_no:   cell number as specified in cell_str, crated in the cell file manager
%monkey:    self-explanatory, enter as string
%area:      idem
%OUTPUTS:
%RT_out:    cell with RTs in each condition
%conditions_plotted_names: name of each condition in RF-centered coordinates
%conditions_plotted_for_computation: name of each condition in absolute coordenates 


base_dir='/Freiwald/ppolosecki/lspace';
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

if ismember(cell_no,bad_files(monkey,area))
    error('This is one of them bad cells, kid. We don''t want no bad cells in here')
end

load(cell_file);
load(results_file)
load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.mat{end}));

[surf_str] = trial_info([],tds);

%Extract reaction times;

trials_used=strcmp({tds{1}.trials.type}','Normal') & strcmp({tds{1}.trials.outcome}','Success');
RTs=[tds{1}.trials.reactime]';

min_index=find(unique([surf_str.phi]')==results_data(cell_no).closest_surf_phi);
phis=[surf_str.phi]';
brts=[surf_str.brt]';
if any(phis==45); 
    shifted=1;
    phis=phis-45;
    brts=brts-45;
else
    shifted=0;
end




conditions_plotted_for_computation=allcomb(circshift(unique(phis),-min_index+1),circshift(unique(phis),-min_index+1));
conditions_plotted_names=allcomb(unique(phis),unique(phis));

RT_out=cell(size(conditions_plotted_names,1),1);
for i=1:size(conditions_plotted_names,1)
    RT_out{i}=RTs(phis==conditions_plotted_for_computation(i,1) & brts==conditions_plotted_for_computation(i,2) & trials_used);
end

if shifted==1
conditions_plotted_for_computation=conditions_plotted_for_computation+45;
end
