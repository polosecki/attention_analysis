function [tmap] = trialselect_rewvals(surf_str, rewval,tmap)

% Select trials by outcomes (reduce the trialmap)
Rews = [surf_str.rew];
rews = Rews(tmap.idx);
idx =  find(rews ==rewval);
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(idx);
