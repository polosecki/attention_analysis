function [timing]= measure_timing_of_effects_fun(monkey,area,use_high_res_data,noise_model)
if nargin<4
    noise_model='poisson';
end

%max_align=3; %set to 3 if BRT aligned stuff exists %1:align to surf onset; 2:align to saccade onset
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
base_dir='/Freiwald/ppolosecki/lspace';

cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

load(cell_file)
good_files=true(length(cell_str),1);
%% Load general files
addpath('..')

tic
if use_high_res_data
    load(fullfile(base_dir,'polo_preliminary','attention_analysis','population_GLM',[area '_' monkey '_' noise_model '_highres_GLMs.mat']))
else
    load(fullfile(base_dir,'polo_preliminary','attention_analysis','population_GLM',[area '_' monkey '_' noise_model '_GLMs.mat']))
end
toc

%% measure timing of attention:

%raw_cell_no=25; %Quincy LIP cell 15: has saccadic a signal independent of attention!!!

make_figures=false;
% if ismember(raw_cell_no,bad_files(monkey,area))
%     error('Cell belongs to bad cells')
% end
%
% cell_no=raw_cell_no-sum(bad_files(monkey,area)<raw_cell_no);

attention_long_enough=nan(length(cell_str),1);
attention_duration=nan(length(cell_str),1);
attention_onset=nan(length(cell_str),1);
saccade_onset=nan(length(cell_str),1);
saccade_long_enough=nan(length(cell_str),1);


for raw_cell_no=1:length(attention_long_enough)
    
    if ~ismember(raw_cell_no,bad_files(monkey,area))
    cell_no=raw_cell_no-sum(bad_files(monkey,area)<raw_cell_no);
        %attn_sig=all_cell_results{cell_no,1}.GLM(1).Fsig(4,:);
        %attn_sig=all_cell_results{cell_no,1}.GLM(2).Fsig(8,:);
        attn_sig=all_cell_results{cell_no,1}.GLM(1).Fsig(1,:);
        t=all_cell_results{cell_no,1}.time;
        [t0,t0_ind]=min(abs(t-0.03));
        t0=t(t0_ind);p_thres=0.05;
        
        if attn_sig(t0_ind)<p_thres
            test_t=t(t<=t0);
            test_sig=attn_sig(t<=t0);
            timing_indx=find(test_sig>p_thres,1,'last');
            timing_indx=timing_indx+1;
            if isempty(timing_indx);
                timing_indx=1;
            end
            %t_split=test_t(timing_indx);
            is_attn_long_enough=true;
            attn_duration=-1;
        elseif attn_sig(t0_ind)>p_thres
            test_t=t(t>t0);
            test_sig=attn_sig(t>t0);
            timing_indx=find(test_sig<p_thres,1,'first');
            diff_time=test_t(2:end);
            if timing_indx~=1
                is_attn_long_enough=diff(diff_time(find(diff(test_sig>p_thres),2,'first')))>0.25;
                attn_duration=diff(diff_time(find(diff(test_sig>p_thres),2,'first')));
                if isempty(attn_duration);attn_duration=0;end
            else
                is_attn_long_enough=diff(diff_time([1 find(diff(test_sig>p_thres),1,'first')]))>0.25;
                attn_duration=diff(diff_time([1 find(diff(test_sig>p_thres),1,'first')]));
            end
            
            
            
            
            
        end
        t_split_attn=test_t(timing_indx)
        
        if isempty(is_attn_long_enough);is_attn_long_enough=true;end
        attention_long_enough(raw_cell_no)=is_attn_long_enough;
        attention_duration(raw_cell_no)=attn_duration;
        attention_onset(raw_cell_no)=t_split_attn;
        
        if make_figures
            semilogy(t,attn_sig,'.-')
            line(xlim,[.05 .05],'Color','black','LineStyle','--')
            line([0 0],ylim,'Color','black')
        end
        
        %load struct files by garch.m
        %/Freiwald/ppolosecki/lspace/polo_preliminary/attention_analysis/population_GLM/garch.m
        
        % Measure timing of saccade tuning:
        
%        sacc_sig=all_cell_results{cell_no,2}.GLM(2).Fsig(5,:);
        sacc_sig=all_cell_results{cell_no,2}.GLM(2).Fsig(2,:);
        t=all_cell_results{cell_no,2}.time;
        [t0,t0_ind]=min(abs(t+0.05));
        t0=t(t0_ind);
        p_thres=0.05;
        
        if sacc_sig(t0_ind)<p_thres
            test_t=t(t<=t0);
            test_sig=sacc_sig(t<=t0);
            timing_indx=find(test_sig>p_thres,1,'last');
            %number of consecutive bins where this happens:
            %length(test_sig)-find(diff(test_sig>p_thres),1,'last')
            if timing_indx~=length(test_t)
                timing_indx=timing_indx+1;
                t_split_sacc=test_t(timing_indx)
            else
                t_split_sacc=t0;
                
            end
            is_sacc_long_enough=abs(t_split_sacc)>=0.1;
            
        elseif sacc_sig(t0_ind)>p_thres
            %     test_t=t(t>t0);
            %     test_sig=sacc_sig(t>t0);
            %     timing_indx=find(test_sig<p_thres,1,'first');
            %     t_split_sacc=test_t(timing_indx)
            t_split_sacc=nan;
            is_sacc_long_enough=false;
        end
        
        saccade_onset(raw_cell_no)=t_split_sacc;
        saccade_long_enough(raw_cell_no)=is_sacc_long_enough;
        
        if make_figures
            figure;
            semilogy(t,sacc_sig,'.-')
            line(xlim,[.05 .05],'Color','black','LineStyle','--')
            line([0 0],ylim,'Color','black')
        end
        
    end
end


timing.monkey=monkey;
timing.area=area;
timing.attention.onset=attention_onset;
timing.attention.duration=attention_duration;
timing.saccade.onset=saccade_onset;