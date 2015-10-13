function [p_val] = test_dicrimination_timing(sacc_dir,brt_dir,ok_trials,n_perms)
% INPUTS
% sacc_dir: vector of saccade directions
% BRT_dir:  vector of BRT direction
%n_perms: number of permutations to consider (determines precision of sig levels)
% OUTPUT
% p_val: fraction of permutations that give the number of trials where RT
% occurs after BRT turns out larger than the observed one

p_val = 0;

brt_dir = brt_dir(ok_trials);
sacc_dir = sacc_dir(ok_trials);
abs_val = unique(brt_dir);

brt_dir(brt_dir==abs_val(2))=abs_val(1);
brt_dir(brt_dir==abs_val(4))=abs_val(3);
sacc_dir(sacc_dir==abs_val(2))=abs_val(1);
sacc_dir(sacc_dir==abs_val(4))=abs_val(3);

