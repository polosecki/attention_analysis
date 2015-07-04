function p = heiko_log_surf_params(td,nocue)

% Fetch surface parameters of the cued surface
p.cuedsurfnames = cell(1,length(td.trials));
p.cuedsurfangles = zeros(1,length(td.trials));
p.cuedbitdurs = cell(1,length(td.trials));
p.cuedbitdirs = cell(1,length(td.trials));
p.ncuedsurfnames = cell(1,length(td.trials));
p.ncuedsurfangles = zeros(1,length(td.trials));
p.ncuedbitdurs = cell(1,length(td.trials));
p.ncuedbitdirs = cell(1,length(td.trials));
p.flickerdur = zeros(1,length(td.trials));
p.pausedur = zeros(1,length(td.trials));
p.bit_total_time = zeros(1,length(td.trials));
for tk=1:length(td.trials)
   if strcmp( td.trials(tk).type, 'Normal') %scs added this to take care of catch trials 3-29-13
    p.cuedsurfnames{tk} = td.trials(tk).targsurf.name;
    if isfield(td.trials(tk),'distsurf')
        p.ncuedsurfnames{tk} = td.trials(tk).distsurf(1).name;
    end
   end %scs
   if strcmp(td.trials(tk).type, 'Catch')%scs
       if td.trials(tk).cueseq == [0 0]
           p.cuedsurfnames{tk} = td.trials(tk).distsurf(1).name;
           p.ncuedsurfnames{tk} = td.trials(tk).distsurf(2).name; 
       end
   end %scs
   
        if strcmp( td.trials(tk).type, 'Normal')%scs
    p.bit_total_time(tk) = 1000*sum(td.trials(tk).targsurf.bitdurs)/td.framerate;
        end %scs
        if strcmp(td.trials(tk).type, 'Catch')
           p.bit_total_time(tk) = nan;
          
        end %scs
        
    % Find the position of the cued surface
     if strcmp( td.trials(tk).type, 'Normal')%scs
    allnames = [td.doc_data(tk).FIELDPARAMS_HEIKO.nameField];
    idx = find(strcmp(allnames,p.cuedsurfnames{tk}));
    p.cuedsurfangles(tk) = str2double(td.doc_data(tk).FIELDPARAMS_HEIKO(idx).posPhi{1});
    if nocue
        p.flickerdur(tk) = 0;
        p.pausedur(tk) = 0;
    else
        p.flickerdur(tk) = str2double(td.doc_data(tk).FIELDPARAMS_HEIKO(idx).flickerDur{1});
        p.pausedur(tk) = str2double(td.doc_data(tk).FIELDPARAMS_HEIKO(idx).pauseDur{1});
    end
    if ~isempty(p.ncuedsurfnames{tk})
        idx = find(strcmp(allnames,p.ncuedsurfnames{tk}));
        p.ncuedsurfangles(tk) = str2double(td.doc_data(tk).FIELDPARAMS_HEIKO(idx).posPhi{1});
    end
     end
    if strcmp( td.trials(tk).type, 'Normal')%scs
    p.cuedbitdurs{tk} = td.trials(tk).targsurf.bitdurs;
    p.cuedbitdirs{tk} = td.trials(tk).targsurf.bitdirs;
    p.ncuedbitdurs{tk} = td.trials(tk).distsurf(1).bitdurs;
    p.ncuedbitdirs{tk} = td.trials(tk).distsurf(1).bitdirs;
    end %scs
end
