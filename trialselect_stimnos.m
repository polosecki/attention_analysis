function [tmap] = trialselect_stimnos(tmap,stimnos)

% Select trials by stimulus numbers (reduce the trialmap)

idx = find(ismember(tmap.sn,stimnos));
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(idx);
