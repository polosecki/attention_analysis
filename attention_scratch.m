%imagesc(grand_psth.to_saccade_onset);

phi=[surf_str.phi];
brt=[surf_str.brt];
[sorted_phi,I_sphi]=sort(phi);
temp=grand_psth.to_stim_onset(I_sphi,:);
imagesc(temp(trials_used(I_sphi),:));
%%

%close all
chosen_brt=135;
good_brt=([surf_str.brt]==chosen_brt)';
%good_phi=(ismember([surf_str.phi],[315]))';
good_phi=(~ismember([surf_str.phi],[135]))';
%good_phi=(~ismember([surf_str.phi],[135 315]))';
%good_phi=(~ismember([surf_str.phi],[45 225]))'
rt=[surf_str.rts];
[sorted_rts,I_srts]=sort(rt);
temp2=grand_psth.to_saccade_onset(I_srts,:);

data_seen=temp2(trials_used(I_srts) & good_brt(I_srts) & good_phi(I_srts),:);

ker=ones(5,6); ker=ker./sum(ker(:));
figure; imagesc(data_seen)
%colormap('gray')

data_seen=conv2(data_seen,ker,'same');
figure;
imagesc(data_seen);
%colormap('gray')

figure; plot(sorted_rts(trials_used(I_srts) & good_brt(I_srts) & good_phi(I_srts)));

figure; plot(nanmean(data_seen(end-15:end,:)));hold on
plot(nanmean(data_seen(1:15,:)),'r')