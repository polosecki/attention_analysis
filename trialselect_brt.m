function [tmap] = trialselect_brt(tmap, selbrt, tds)
%scs 12/06/28
% function [tmap] = trialselect_targets(tmap,target,t_idx)
% %scs 12/06/28

% Select trials by outcomes (reduce the trialmap)
%find the brt's
for j= 1:length(tds{1}.trials)
    if ~strcmp(tds{1}.trials(j).type,'Normal') brt(j) = nan; 
    else
    brt(j)= tds{1}.trials(j).targsurf.brtdir;
    end 
end

idx = find(ismember(brt,selbrt));
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(ismember(tmap.idx,idx)); 
