function [surf_str, surf_str_extra, td] = get_behavioral_data(monkey,area,cell_no)
%monkey='Michel';%'Quincy';%
%area='PITd';
%cell_no=15;


base_dir=fullfile('/Freiwald/ppolosecki/lspace',lower(monkey));

cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);

load(fullfile(cell_file));

load(fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.mat{end}));
ifname=fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.hdf{end});

% ------------ END OF USER-EDITABLE BLOCK -----------------
addpath('..')



fid = dhfun('open',ifname);
tmap = dh_get_trialmap_struct(fid);
dhfun('close',fid);
%% Makes all relevant data structures
%trls = [tds{1}.trials]; %holds all the basic trial parameters
%this gets you the name, number, rho and phi for the target surface of
%every trial
%[surf_str] = surface_info_withcatch(tmap,tds);
[surf_str] = trial_info(tmap,tds);
surface_numbers = [surf_str.numb];

%This is Michael Code (partially redundant):----

td = tds{1};
nocue = false;
surf_str_extra = heiko_log_surf_params(td,nocue); % used to be called p in michael code
%---
%Sanity checks:
if any(strcmp({surf_str.out},'Success')'~=(tmap.oc==7))
    error('Tmap and surface info are not consistent')
end
test1=surf_str_extra.cuedsurfnames;
test2=[surf_str.name];
if ~all(strcmp(test1,test2))
    error('The surf_str surface info structures are not consistent')
end
%-----
