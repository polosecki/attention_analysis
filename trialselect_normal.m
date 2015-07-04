function [tmap] = trialselect_normal(normal,tmap)

% Select trials by outcomes (reduce the trialmap)

idx =  find(normal(tmap.idx) == 1);
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(idx);
