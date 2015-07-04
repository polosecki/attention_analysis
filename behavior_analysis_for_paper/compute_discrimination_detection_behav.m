function [results] = compute_discrimination_detection_behav(monkey,area,cell_no)
%monkey='Quincy';%'Michel';%
%area='PITd';
%cell_no=29;

[surf_str, surf_str_extra, td] = get_behavioral_data(monkey,area,cell_no);

outcomes={surf_str.out}';
RTs=[surf_str.rts]';
normal_tr=strcmp({surf_str.type}','Normal');

allowed_outcomes = {'Early reaction', 'Wrong target', 'Wrong surface target',...
    'Success', 'No target response'};

used_trials=normal_tr & ismember(outcomes,allowed_outcomes);
nr = used_trials & normal_tr & ismember(outcomes,'No target response');


a=td.doc_data;
b=[a.TARGETSPARAMS];
if any(~strcmp(unique([b.actDelay]),'0')) || unique(surf_str_extra.pausedur(used_trials))~=1 ...
        || unique(surf_str_extra.flickerdur)~=0
    error('RT not measured relative to BRT')
end
clear a b

frame_dur=1/td.framerate; % in seconds
BRT_time=cellfun(@sum,surf_str_extra.cuedbitdurs)'*(1/td.framerate);
nBRT_time=cellfun(@sum,surf_str_extra.ncuedbitdurs)'*(1/td.framerate);
%RTs=RTs+BRT_time;

BRT_time(~used_trials)=nan;
nBRT_time(~used_trials)=nan';
RTs(~used_trials)=nan';
%keyboard
td.visiko_versions
sacc_dir=[td.trials.actsacangle]';
sacc_dir(~used_trials)=nan;
brt_dir=[surf_str.brt]';
brt_dir(~used_trials)=nan;

nbrt_dir=nan(size(brt_dir));
a=td.trials;
b=[a(used_trials).distsurf];
nbrt_dir(used_trials)=[b.brtdir]';
clear a b

%RT measured with respect to trial onset:
k = unique(surf_str_extra.pausedur(used_trials));
RT_absolute=RTs + k + BRT_time;
RT_rel_dist = RT_absolute - k - nBRT_time;
%Calculateting RTs from hdf:
% fid = dhfun('open',ifname);
% out_RTs = calculate_RTs_from_hdf(fid,tmap);
% dhfun('close',fid);

%Calculating detection chance levels:
early_tr=used_trials & ismember(outcomes,{'Early reaction'});
%p_val_detection = binomialcdf(sum(early_tr),sum(used_trials),0.5);
%p_val_detection2 = binomialcdf(sum(RT_rel_dist(used_trials)<0),sum(used_trials),0.5);
n_perms=100000;
%min(RT_absolute(used_trials))
RT_absolute(nr)=0;
[p_val_detection_target, targ_detect_chance_level, targ_detect_conf_int] = test_discrimination_timing(RT_absolute(used_trials),BRT_time(used_trials)+k,n_perms)
if p_val_detection_target==0; p_val_detection_target = 1/ n_perms; end
percent_target_detected = sum(RTs(used_trials)>0 & RTs(used_trials)<2.7)/sum(used_trials) * 100



[p_val_detection_distractor, dist_detect_chance_level, dist_detect_conf_int] = test_discrimination_timing(RT_absolute(used_trials),nBRT_time(used_trials)+k,n_perms)
if p_val_detection_distractor==0; p_val_detection_distractor = 1/ n_perms; end
percent_dist_detected  = sum(RT_rel_dist(used_trials)>0 & RT_rel_dist(used_trials)<2.7)/sum(used_trials) * 100
%keyboard

%Calculating discrimination chance levels:
RTs(nr)=-0.1;
target_discrim_behavior = (RTs > 0 & RTs < 2.7) & used_trials & ~isnan(sacc_dir);
ok_tr1=target_discrim_behavior & (brt_dir ~= sacc_dir); %ismember(outcomes,{'Wrong target', 'Wrong surface target'});
p_val_discrimination_target = binomialcdf(sum(ok_tr1),sum(used_trials),0.5)
percent_target_discriminated = (1-sum(ok_tr1)/sum(used_trials)) * 100
target_discrim_chance_conf_int = binoinv([0.05 0.95],sum(used_trials),0.5)/sum(used_trials) * 100

RT_rel_dist(nr)=-0.1;
trials_distractor_discrim= (RT_rel_dist > 0 & RT_rel_dist< 2.7) & used_trials & ~isnan(sacc_dir);
distractor_discrim_behavior = (nbrt_dir ~= sacc_dir) & trials_distractor_discrim;
p_val_discrimination_dist = binomialcdf(sum(distractor_discrim_behavior),sum(trials_distractor_discrim),0.5)
percent_distractor_discriminated = (1-sum(distractor_discrim_behavior)/sum(trials_distractor_discrim)) * 100

results.performance.detection_target=percent_target_detected;
results.performance.detection_distractor=percent_dist_detected;
results.performance.discrimination_target=percent_target_discriminated;
results.performance.discrimination_distractor=percent_distractor_discriminated;

results.chance.detection_target_median = targ_detect_chance_level;
results.chance.detection_target_confidence = targ_detect_conf_int;
results.chance.detection_dist_median = dist_detect_chance_level;
results.chance.detection_dist_confidence = dist_detect_conf_int;
results.chance.discrimination_confidence = target_discrim_chance_conf_int;
results.distributions.RT_abs = RT_absolute(used_trials);
results.distributions.BRT_abs = BRT_time(used_trials)+k;
results.distributions.nBRT_abs = nBRT_time(used_trials)+k;

results.sig.detection_target=p_val_detection_target;
results.sig.detection_distractor=p_val_detection_distractor;
results.sig.discrimination_target=p_val_discrimination_target;
results.sig.discrimination_distractor=p_val_discrimination_dist;
