function [ts,te,inclidx] = find_ref_points_heiko_PP(fid,tmap,start_trig,end_trig)
% Find the intervals of interest based on trigger points and offsets
%
%Adapted from Michael Borisov's find_ref_points.m by Pablo Polosecki (8/1/2014)
%
%INPUTS:
%fid: the file ID
%tmap: tmap provided by the function dh_get_trialmap_struct
%start_trig: struct with fields 'type', 'name' and 'offset'
%end_trig: analogous, but for the end trigger
%surf_str: as provided by the trial_info function
%
%OUTPUTS:
%ts: vector of start times (in ns)
%te: vector of end times (in ns)
%inclidx: vector with is true for trials with satisfactorily calculated reference
%points


global DH 

excl = zeros(length(tmap.ts),1);

if strcmp(start_trig.type, 'TrialStart')
    % Just use the trialstart time which is present for all trials
    ts = tmap.ts;

elseif strcmp(start_trig.type,'TrialEnd')
    % Just use the trialend time which is present for all trials
    ts = tmap.te;

elseif strcmp(start_trig.type,'IntervalBegin')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,start_trig.name);
    [ts,excl] = proc_interval(tmap,its,ite,its,excl);
elseif strcmp(start_trig.type, 'IntervalEnd')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,start_trig.name);
    [ts,excl] = proc_interval(tmap,its,ite,ite,excl);
elseif strcmp(start_trig.type,'Marker')
    % Find markers which are contained by trials
    mt = dhfun('GETMARKER',fid,start_trig.name);
    if strcmp(start_trig.name,'Fixation_out') %  Makes up for multiple fix out confusion
        stim_on_times = dhfun('GETMARKER',fid,'Stimulus_onset');
        [ts,excl] = proc_marker_fix_out(tmap,mt,stim_on_times,excl);
    else
        [ts,excl] = proc_marker(tmap,mt,excl);
    end
else
    error(['Unknown start_trig.type specification' start_trig.type]);
end

if strcmp(end_trig.type, 'TrialStart')
    te = tmap.ts;
    
elseif strcmp(end_trig.type, 'TrialEnd')
    te = tmap.te;
    
elseif strcmp(end_trig.type,'IntervalBegin')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,end_trig.name);
    [te,excl] = proc_interval(tmap,its,ite,its,excl);
elseif strcmp(end_trig.type,'IntervalEnd')
    % Find intervals which are at least partly contained by trials
    [its,ite] = dhfun('GETINTERVAL',fid,end_trig.name);
    [te,excl] = proc_interval(tmap,its,ite,ite,excl);
elseif strcmp(end_trig.type,'Marker')
    % Find markers which are contained by trials
    mt = dhfun('GETMARKER',fid,end_trig.name);
    if strcmp(end_trig.name,'Fixation_out') %  Makes up for multiple fix out confusion
        stim_on_times = dhfun('GETMARKER',fid,'Stimulus_onset');
        [te,excl] = proc_marker_fix_out(tmap,mt,stim_on_times,excl);
    else
        [te,excl] = proc_marker(tmap,mt,excl);
    end
   
else
    error(['Unknown end_trig.type specification: ' end_trig.type]);
end

% Apply the offsets
ts = ts + round(1e9*start_trig.offset);
te = te + round(1e9*end_trig.offset);

% Exclude also the trials where te < ts
excl(find(~excl & (te < ts))) = 1;

% Now, exclude all trial intervals which were marked for exclusion
if any(excl)
    inclidx = find(~excl);
    ts = ts(inclidx);
    te = te(inclidx);
else
    inclidx = 1:length(ts);
end
% End of main function
% --------------------

function [tp,excl] = proc_interval(tmap,its,ite,itp,excl)
tp = zeros(length(tmap.ts),1);
for i=1:length(tmap.ts)
    idx = find(ite >= tmap.ts(i) & its <= tmap.te(i));
    % If there are no intervals found within this trial, or
    % more than one interval, exclude this trial
    if length(idx) ~= 1
        excl(i) = 1;
    else
        tp(i) = itp(idx);
    end
end

function [tp,excl] = proc_marker(tmap,mt,excl)
tp = zeros(length(tmap.ts),1);
for i=1:length(tmap.ts)
    idx = find(tmap.ts(i) <= mt & tmap.te(i) >= mt);
    % If there are no markers found within the trial, or
    % more than one marker, exclude this trial
    if length(idx) ~= 1
        excl(i) = 1;
    else
        tp(i) = mt(idx);
    end
end


function [tp,excl] = proc_marker_fix_out(tmap,fix_out_times,stim_on_times,excl)
tp = zeros(length(tmap.ts),1);
for i=1:length(tmap.ts)
    %Look for the stimulus onset time
    stim_on_id=find(tmap.ts(i) <= stim_on_times & tmap.te(i) >= stim_on_times);
    if length(stim_on_id) ~=1
        excl(i) = 1; %if there is more than one stim onset, exclude trial
    else
        this_stim_on_time=stim_on_times(stim_on_id);
        idx = find(this_stim_on_time <= fix_out_times & tmap.te(i) >= fix_out_times);%look for fix out events btween stim onset and trial end
        if length(idx) < 1 %if there is none, exlcude
            excl(i) = 1;
        elseif length(idx)>1 %if there are many, take the first one
            disp(['Fixed trial ' num2str(i)])
            tp(i)=fix_out_times(idx(1));
        else
            tp(i) = fix_out_times(idx);%if there is one then :)
        end
    end
end