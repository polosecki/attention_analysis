function [t_brt,inclidx] = find_BRT_times_heiko_PP(fid,tmap,use_distractor)

trial_start=tmap.ts;
trial_end=tmap.te;

t_brt=nan(length(trial_start),1);
inclidx=true(size(t_brt));

if nargin>2 && use_distractor
all_brt_times = dhfun('GETMARKER',fid,'BRT_target_onset');
else
all_brt_times = dhfun('GETMARKER',fid,'BRT_target_onset');
end

for i=1:length(t_brt)
    stim_on_id=find(tmap.ts(i) <= stim_on_times & tmap.te(i) >= stim_on_times);
    this_trial_brt_indx=find(trial_start(i)<=all_brt_times & trial_end(i) >= all_brt_times);
    if length(this_trial_brt_indx) > 1
        error(['More than one BRT onset found in trial ' num2srt(i)])
    elseif length(this_trial_brt_indx)==1
        t_brt(i)=all_brt_times(this_trial_brt_indx);
    else
        inclidx(i)=false;
    end
end
    
