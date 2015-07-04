function [surf_str] = surface_info_notmap(tds)%surface_info_new(tmap, tds)
% to get the configuration of stimuli for a given trial on the attention
% task
%this doesn't require the tmap, so it can be used just for behavior files.

surf_names = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.nameField];
surf_phis= [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.posPhi];
surf_rhos = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.posRho];
num_phi = 1:length(unique(surf_phis))
uni_phi = str2double(unique(surf_phis))
coh_lev = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.brtCoher];
trls = tds{1}.trials;

for m = 1:length(trls); 
   %surf_str(m).idx =  tmap(1).idx(m); %put back in!!!! 3/15/13
    if strcmp(tds{1}.trials(m).type,'Neutral')
        surf_str(m).name = nan;
        surf_str(m).numbtmap = nan;
        surf_str(m).numb = nan;
        surf_str(m).rho = nan;
        surf_str(m).brt = nan;
        surf_str(m).phi = nan;
        surf_str(m).coh = nan; 
        surf_str(m).rew = nan;
        surf_str(m).type = ('Neutral');
    end %need to remove neutral phi to keep from getting a surface 3 
    if strcmp(tds{1}.trials(m).type, 'Catch') 
            if tds{1}.trials(m).cueseq == [0 0];
    cued_surf =  tds{1}.trials(m).distsurf.name; 
            end
    surf_str(m).name = cued_surf;
    surf_str(m).rts = nan;
    surf_str(m).type = ('Catch');
    end
   if strcmp(tds{1}.trials(m).type, 'Normal')        
surfname = [];
surfnumb = [];
surfname = [tds{1}.trials(m).targsurf.name];
surf_str(m).brt = tds{1}.trials(m).targsurf.brtdir;
surf_str(m).rts = tds{1}.trials(m).reactime;
surf_str(m).name = surf_names(strcmp(surfname,surf_names)==1);
surf_str(m).type = ('Normal'); 
    end
%surfnumb = tmap(1).sn(m); %put back in, took out for testing SCS 3/15/13
surfnumb = nan;
surfname = [];
surfname = surf_str(m).name;
%surf_str(m).numbtmap =surfnumb; 

surf_str(m).phi = str2double(surf_phis(strcmp(surfname,surf_names)==1));
%%commented 2/6/13 scs, something funny, check
%FIX THIS
%if num_phi>2
%surf_str(m).numb = num_phi(strcmp(surfname, surf_names)==1); %fix to
%number based on phi %this might be fixed in this new version (130219 scs)
surf_str(m).numb = num_phi(surf_str(m).phi== uni_phi);
surf_str(m).numbtmap =surfnumb;
% else
%     surf_str(m).numb =surfnumb; 
% end
% surf_str(m).numb = num_phi(strcmp(surf_phis)
surf_str(m).rho = str2double(surf_rhos(strcmp(surfname,surf_names)==1));

surf_str(m).coh = str2double(coh_lev(strcmp(surfname,surf_names)==1));
surf_str(m).rew = tds{1}.trials(m).reward;
surf_str(m).out = tds{1}.trials(m).outcome;


    end
