clear all
monkey = 'Michel';
area = 'LIP';
addpath('/Freiwald/ppolosecki/lspace/polo_preliminary/attention_analysis/')
included_sessions={};
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
load(cell_file)

basedir=fullfile('/Freiwald/ppolosecki/lspace', lower(monkey));
[curr_dir, ~, ~]=fileparts(mfilename('fullpath'));


for cell_no=3:3%2:2%4:11 %2
    if ~ismember(cell_no,bad_files(monkey,area))
        cell_no
        log_dir=fullfile(basedir,cell_str(cell_no).dir);
        proc_dir= fullfile(log_dir,'proc');
        matfname=fullfile(proc_dir,cell_str(cell_no).attention.mat{1});
        [~,fbody,~]=fileparts(cell_str(cell_no).attention.mat{1});
        log_file=fullfile(log_dir,[fbody '.log']);
        
        cd('/Freiwald/ppolosecki/lspace/BorisovMatlab/proc_visiko_lrsac/visiko_log_parser')
        [tds,tsls,status] = visiko_log_parser_fun(log_file,0,10000,true);
        %[task_datas,task_start_lines,status] = visiko_log_parser_fun(input_filename,sessions_to_skip,max_sessions_to_read,handle_incomplete_files)
        if ~isempty(status)
            disp('FAIL');
            fprintf('Parse error was detected, possibly incomplete file\n%s\n',...
                status.message);
        else
            disp('OK');
            if ~exist([matfname '.old'],'file')
                movefile(matfname,[matfname '.old'])
            %else
            %    warning('Old backup already exists')
            end
            save(matfname,'tds','tsls');
        end
    end
end
cd(curr_dir)