ifname = '/Freiwald/ppolosecki/lspace/quincy/130509Quincy/proc/Quincy_130509_0041.mat';
%-----------------------------------------

seloutcomes = {'Success' 'Wrong target' 'Wrong surface target' 'Early reaction'};
seltrialtypes = {'Normal'};
% Load the log/mat file
if ~exist('td','var')
    tds = {};
    tsls = '';
    load(ifname);
    td = tds{1};
    clear tds
    clear tsls
end

outcomes = {td.trials.outcome};
trialtypes = {td.trials.type};
reactimes = [td.trials.reactime];

cond_outcomes = ismember(outcomes,seloutcomes);
cond_trialtypes = ismember(trialtypes,seltrialtypes);
cond = cond_outcomes & cond_trialtypes;
idx = find(cond);

framerate = td.framerate;

% Compute the attended surface BRT time based on frame log
react = zeros(1,length(idx));
adjreact = zeros(1,length(idx));
tbits_vect = zeros(1,length(idx));
for i=1:length(idx)
    ta = td.trials(idx(i));
    tbits = sum(ta.targsurf.bitdurs)/framerate;
    targ_surf_name = ta.targsurf.name;
    all_surf_params = td.doc_data(idx(i)).FIELDPARAMS_HEIKO;
    targ_surf_idx = find(strcmp([all_surf_params.nameField],targ_surf_name),1);
    targ_surf_params = all_surf_params(targ_surf_idx);
    flickerdur = sscanf(targ_surf_params.flickerDur{1},'%f');
    if flickerdur~=0
        disp(['flickerdur = ' mun2str(flickerdur)])
        disp(['i= ' num2str(i) ', idx=' num2str(idx(i))])
    end
    pausedur = sscanf(targ_surf_params.pauseDur{1},'%f');
    if pausedur~=1
        disp(['pausedur = ' mun2str(pausedur)])
        disp(['i= ' num2str(i) ', idx=' num2str(idx(i))])
    end
    tbrt = flickerdur+pausedur+tbits;
    tbits_vect(i) = tbits;
    % Check that the reaction time + tbrt >= 0
    react(i) = reactimes(idx(i));
    adjreact(i) = tbrt + reactimes(idx(i));
end
figure(1);
plot(react);
title('Reaction times relative to target BRT');
figure(2);
plot(adjreact);
title('Reaction times relative to stimulus onset');

if(any(adjreact < 0))
    error('Bad trial(s) found');
end
