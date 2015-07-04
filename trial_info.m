function [surf_str] = trial_info(tmap,tds)%surface_info_new(tmap, tds)

%tmap argument is optional

% to get the configuration of stimuli for a given trial on the attention
% task
%this doesn't require the tmap, so it can be used just for behavior files.

surf_names = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.nameField];
surf_phis= [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.posPhi];
surf_rhos = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.posRho];
num_phi = 1:length(unique(surf_phis));
uni_phi = str2double(unique(surf_phis))
coh_lev = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.brtCoher];
trls = tds{1}.trials;

for m = 1:length(trls);
    %surf_str(m).idx =  tmap(1).idx(m); %put back in!!!! 3/15/13
    surf_str(m).type = tds{1}.trials(m).type; % 'Catch', 'Normal' or 'Neutral'
    if ~strcmp(tds{1}.trials(m).type, 'Normal') %Non-normal trials
        if tds{1}.trials(m).cueseq == [0 0];
            cued_surf =  tds{1}.trials(m).distsurf.name;
        else
            error('Not ready to handle cue-switch trials')
        end
        surf_str(m).name = cued_surf;
        surf_str(m).rts = nan;
        surf_str(m).brt = tds{1}.trials(m).distsurf.brtdir;
            
    else %normal trials
        surfname = [tds{1}.trials(m).targsurf.name];
        surf_str(m).brt = tds{1}.trials(m).targsurf.brtdir;
        surf_str(m).rts = tds{1}.trials(m).reactime;
        surf_str(m).name = surf_names(strcmp(surfname,surf_names)==1);
    end
    if ~isempty(tmap)
        surfnumb = tmap.sn(m);
    end
    surfname = surf_str(m).name;
    
    surf_str(m).phi = str2double(surf_phis(strcmp(surfname,surf_names)));
    surf_str(m).numb = num_phi(surf_str(m).phi== uni_phi);
    if ~isempty(tmap)
    surf_str(m).numbtmap =surfnumb;
    end
    surf_str(m).rho = str2double(surf_rhos(strcmp(surfname,surf_names)==1));
    
    surf_str(m).coh = str2double(coh_lev(strcmp(surfname,surf_names)==1));
    surf_str(m).rew = tds{1}.trials(m).reward;
    surf_str(m).out = tds{1}.trials(m).outcome;
    
    
end
