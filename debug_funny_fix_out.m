trials_used; % Succesful, normal trials with good markers
ideal_trials=strcmp({surf_str.out},'Success')' & strcmp({surf_str.type},'Normal')';%succesful normal trials
strange_trials=find(ideal_trials & ~trials_used);%succesful normal trials with funny markers
typical_trials=find(trials_used);%same as trials_sued, but index-based


fid = dhfun('open',ifname,'r'); %open hdf file

dhfun('getfidinfo',fid)
%time start of first strange trial

period=dhfun('getcontsampleperiod',fid,129);%sample period of x-eyepositiont
%[nsamp,nchan]=dhfun('getcontsize',fid,129);


%%
tr_i=11; % get a trial index
sample_start=floor(tmap.ts(strange_trials(tr_i))/period);
sample_end=floor(tmap.te(strange_trials(tr_i))/period);

use_good_trials=1;
if use_good_trials
    sample_start=floor(tmap.ts(typical_trials(tr_i))/period);
    sample_end=floor(tmap.te(typical_trials(tr_i))/period);
else
    sample_start=floor(tmap.ts(strange_trials(tr_i))/period);
    sample_end=floor(tmap.te(strange_trials(tr_i))/period);
end
output1 = dhfun('readcont',fid,129,sample_start,sample_end,1,1);
output2 = dhfun('readcont',fid,130,sample_start,sample_end,1,1)
%Time axes for plotting, has bad length, why??
%t_axis=0:period:(length(output1))*period;
plot(output1); hold on
plot(output2,'r')
%%
dhfun('close',fid);