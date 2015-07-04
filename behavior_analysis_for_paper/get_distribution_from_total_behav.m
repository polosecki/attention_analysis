function [out_RT_abs, out_BRT_abs, out_nRT_abs, out_nBRT_abs] = get_distribution_from_total_behav(total_behav,used_lim,monkey_used,invert_selection)

if nargin<4
    invert_selection=false;
end

monkeys={'Quincy','Michel'};
if monkey_used~=3
    idx_used = strcmp({total_behav.monkey}',monkeys{monkey_used});
else
    idx_used = true(length(total_behav),1);
end
c={total_behav.distributions};
RT_abs_tot=[];
BRT_abs_tot=[];
nBRT_abs_tot=[];

for i=1:length(c)
    if idx_used(i)
        RT_abs_tot = vertcat(RT_abs_tot,c{i}.RT_abs-1);
        BRT_abs_tot = vertcat(BRT_abs_tot,c{i}.BRT_abs-1);
        nBRT_abs_tot = vertcat(nBRT_abs_tot,c{i}.nBRT_abs-1);
    end
end

good_idx=RT_abs_tot >0;

if invert_selection==false
    high_lim(1)=prctile(BRT_abs_tot,90);
    low_lim(1)=prctile(BRT_abs_tot,10);
    high_lim(2)=prctile(nBRT_abs_tot,90);
    low_lim(2)=prctile(nBRT_abs_tot,10);
    switch used_lim
        case 'low'
            good_idx_BRT = good_idx & BRT_abs_tot<low_lim(1);
            good_idx_nBRT = good_idx & nBRT_abs_tot<low_lim(2);
        case 'high'
            good_idx_BRT = good_idx & BRT_abs_tot>high_lim(1);
            good_idx_nBRT = good_idx & nBRT_abs_tot>high_lim(2);
        otherwise
            good_idx_BRT = good_idx;
            good_idx_nBRT = good_idx;
    end
else
    high_lim=prctile(RT_abs_tot,90);
    low_lim=prctile(RT_abs_tot,10);
    low_lim2=prctile(RT_abs_tot,20);
    switch used_lim
        case 'low'
            good_idx_BRT = good_idx & RT_abs_tot>low_lim & RT_abs_tot<low_lim2;
            good_idx_nBRT = good_idx & RT_abs_tot>low_lim & RT_abs_tot<low_lim2;
        case 'high'
            good_idx_BRT = good_idx & RT_abs_tot>high_lim;
            good_idx_nBRT = good_idx & RT_abs_tot>high_lim;
        otherwise
            good_idx_BRT = good_idx;
            good_idx_nBRT = good_idx;
    end
end

out_RT_abs = RT_abs_tot(good_idx_BRT);
out_BRT_abs = BRT_abs_tot(good_idx_BRT);
out_nRT_abs = RT_abs_tot(good_idx_nBRT);
out_nBRT_abs = nBRT_abs_tot(good_idx_nBRT);