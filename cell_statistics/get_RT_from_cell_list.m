function [RT_data]=get_RT_from_cell_list(monkey,area)
%RT_data is a structure with fields:
%RTs: cell with RTs of correct trials in each condition;
%conditions_RF_centered: name of each condition: column 1 is attended
%surface phi and column 2 is brt direction. In RF-aligned angles.
%conditions_abs_positions: The same, but angles are polar angles in
%absolute space.
%session_index: a session number, that is the same for cells from the same
%day.

cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
load(cell_file)

RT_data=struct;
included_sessions={};
for cell_no=1:length(cell_str)
    if ~ismember(cell_no,bad_files(monkey,area))
        [RT_out, conditions_plotted_names, conditions_plotted_for_computation]=measure_RTs_attention_for_cell(cell_no,monkey,area);
        RT_data(cell_no).RTs=RT_out;
        RT_data(cell_no).conditions_RF_centered=conditions_plotted_names;
        RT_data(cell_no).conditions_abs_positions=conditions_plotted_for_computation;
        if ~ismember(cell_str(cell_no).dir,included_sessions)
            included_sessions=[included_sessions;cell_str(cell_no).dir];
        end
        current_session_index=find(strcmp(cell_str(cell_no).dir,included_sessions));
        RT_data(cell_no).session_index=current_session_index;
    end
    
end

