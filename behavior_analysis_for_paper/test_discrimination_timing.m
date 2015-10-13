function [p_val, chance_level, chance_confidence_interval] = test_discrimination_timing(RT_abs,BRT_abs,n_perms)
% INPUTS
% RT_abs: vector of RTs relative to trials start
% BRT_abs:  vector of BRT onset times relative to trial start
%n_perms: number of permutations to consider (determines precision of sig levels)
% OUTPUT
% p_val: fraction of permutations that give the number of trials where RT
% occurs after BRT turns out larger than the observed one


detected_trials=sum((RT_abs-BRT_abs)>0 & (RT_abs-BRT_abs)<2.7);
rand_perf = nan(n_perms,1);

for i=1:n_perms
    rand_BRT_list = BRT_abs(randperm(length(BRT_abs)));
    rand_perf(i) =sum((RT_abs-rand_BRT_list)>0 & (RT_abs-rand_BRT_list)<2.7);
end

p_val = sum(rand_perf>detected_trials) / length(rand_perf);
chance_level = median(rand_perf)/ length(RT_abs) * 100;
chance_confidence_interval = [prctile(rand_perf,5) prctile(rand_perf,95)] / length(RT_abs) * 100; 
