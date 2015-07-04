function [task_fig] = task_config_notmap(tds) 
%this gets you the name, number, rho and phi for the target surface of
%every trial
% if isempty(ifname);
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\Michel_121116_0028.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\121116_michel_011_r.dh5';
% end
% 
% fid = dhfun('open',ifname);
% tmap = dh_get_trialmap_struct(fid);
trls = [tds{1}.trials]; %holds all the basic trial parameters

%[surf_str] = surface_info_new(tmap, tds); %changed feb 19 2013 scs (from old to new)
[surf_str] = surface_info_withcatch(tds); %for if you have catch trials or no tmap

targ_color = 1; 

num = [surf_str.numb]';
phi = [surf_str.phi]';
rho = [surf_str.rho]';
tgt = [surf_str.brt]';
tmap_sn = [surf_str.numb]';

%nums = 1:length(find(~isnan(unique(phi))));
nums =unique(phi(~isnan(phi)))
if length(nums) >=2
    
for h = 1:length(nums)
    phi_plot(h) = unique(phi(phi==nums(h)));
    rho_plot(h) = max(unique(rho(phi==nums(h))));
    num_plot(h) = nums(h);
    phi_ID(h) = unique(tmap_sn(phi_plot(h)==phi));
end
pp = [phi_ID',phi_plot'];
phi_ind = sortrows(pp, 1)
end


%targets
%rho
targ_rho = str2double(tds{1}.doc_data(1).TARGETSPARAMS.posRho);
targ_phi = unique(tgt(~isnan(tgt))); 

phi_deg = deg2rad(phi_ind(:,2));% for the surfaces
tgt_deg = deg2rad(targ_phi); %for the targets (brt directions)

targ_rhos = [];
figure

targ_rhos(1:length(tgt_deg)) = deal(targ_rho); 

if length(unique(tgt_deg))==8 && targ_color ==1; 
    polar(tgt_deg(1),targ_rhos(1), 'kd')
    hold on
    polar(tgt_deg(2), targ_rhos(2),'c*')
    hold on
    polar(tgt_deg(3), targ_rhos(3),'b*')
    hold on
    polar(tgt_deg(4), targ_rhos(4),'k*')
    hold on
    polar(tgt_deg(5), targ_rhos(5),'g*')
    hold on
    polar(tgt_deg(6), targ_rhos(6),'y*')
    hold on
    polar(tgt_deg(7), targ_rhos(7),'r*')
    hold on
    polar(tgt_deg(8), targ_rhos(8),'m*')
    hold on
else

    polar(tgt_deg, targ_rhos', 'k*')
hold on
end
mark_size=140;
h=polar(phi_deg(1), rho_plot(1),'b.')
set(h,'MarkerSize',mark_size)
hold on
h=polar(phi_deg(2), rho_plot(2), 'r.')
set(h,'MarkerSize',mark_size)
hold on
if length(nums)>2
    h=polar(phi_deg(3), rho_plot(3), 'g.')
    set(h,'MarkerSize',mark_size)
    hold on
    h=polar(phi_deg(4),rho_plot(4), 'm.')
    set(h,'MarkerSize',mark_size)
    hold on
    title('blue is surf 1, red is surf 2, green is 3, magenta is 4')
else
title('blue is surf 1, red is surf 2')
end
hold on




