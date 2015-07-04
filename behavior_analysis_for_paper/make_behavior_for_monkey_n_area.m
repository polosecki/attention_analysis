clear all

monkeys = {'Quincy','Michel'};
areas = {'PITd','LIP'};

for mm=1:2
    for aa=1:2
        monkey = monkeys{mm};%'Michel';
        area = areas{aa};%'LIP';
        
        addpath('..')
        included_sessions={};
        cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
        cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
        load(cell_file)
        
        clear behavior_performance
        for cell_no=1:length(cell_str)
            if ~ismember(cell_no,bad_files(monkey,area))
                if ~ismember(cell_str(cell_no).dir,included_sessions)
                    included_sessions=[included_sessions;cell_str(cell_no).dir];
                    current_session_index=find(strcmp(cell_str(cell_no).dir,included_sessions));
                    behavior_performance(current_session_index) = ...
                        compute_discrimination_detection_behav(monkey,area,cell_no);
                end
            end
        end
        
        save([monkey '_' area '_performance'],'behavior_performance')
    end
end