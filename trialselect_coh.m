function [tmap] = trialselect_coh(tmap, selcoh, tds)
%scs 12/06/28
% function [tmap] = trialselect_targets(tmap,target,t_idx)
% %scs 12/06/28

% Select trials by outcomes (reduce the trialmap)
%find the brt's
% for j= 1:length(tds{1}.trials)
%     if strcmp(tds{1}.trials(j).type,'Neutral') coh(j) = nan; 
%     else
%     brt(j)= tds{1}.trials(j).targsurf.brtdir;
%     end 
% end
for j = 1:length(tds{1}.trials);
if strcmp(tds{1}.trials(j).type,'Neutral') coh(j) = nan; 
    else
surf_name = tds{1}.trials(j).targsurf.name;
all_surfs =[tds{1}.doc_data(1).BLOCKPARAMS_HEIKO_TRIAL.attnSurf];
surf_ind(j) = find(strcmp(surf_name, all_surfs)==1); 
all_cohs = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.brtCoher];
Coh(j) = all_cohs(surf_ind(j));
end
coh= str2double(Coh); 
end

selcohs= find(ismember(coh,selcoh));

idx = find(ismember(tmap.idx,selcohs));
tmap.tn = tmap.tn(idx);
tmap.sn = tmap.sn(idx);
tmap.oc = tmap.oc(idx);
tmap.ts = tmap.ts(idx);
tmap.te = tmap.te(idx);
tmap.idx = tmap.idx(idx);
