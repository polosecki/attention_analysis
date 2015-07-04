function [tmap] = trialselect_outcomes(tmap,outcomes)

% Select trials by outcomes (reduce the trialmap)

idx = find(ismember(tmap.oc,outcomes));
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(idx);
