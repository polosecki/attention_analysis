function [grand_psth_matrix, trials_used, time_axis] = make_grand_psth_matrix(grand_psth_indredients,surf_str,surf_str_extra,td)
%Makes a matrix with individual trial PSTHs:
% Required fields are:
% fid
% tmap
% start_trig
% end_trig
% surf_str
%align %=1 for stim onset or 2 for saccade onset
%spks
%tres
%psthfilter

fid=grand_psth_indredients.fid;
tmap=grand_psth_indredients.tmap;
start_trig=grand_psth_indredients.start_trig;
end_trig=grand_psth_indredients.end_trig;
surf_str=grand_psth_indredients.surf_str;
spks=grand_psth_indredients.spks;
tres=grand_psth_indredients.tres;
psthfilter=grand_psth_indredients.psthfilter;
psthdec=grand_psth_indredients.psthdec;
align=grand_psth_indredients.align;

%Get time start and end (in nanoseconds) for each trial:
%[time_start,time_end,inclidx] = find_ref_points_test(fid,tmap,start_trig,end_trig,surf_str);
%[time_start,time_end,inclidx] = find_ref_points(fid,tmap,start_trig,end_trig);

[time_start,time_end,inclidx] = find_ref_points_heiko_PP(fid,tmap,start_trig,end_trig);

%Prevent future errors!
if find(time_start<0)
    disp('Hey, your time_start is less than zero!')
end

trials_with_good_markers=zeros(length(surf_str),1);
trials_with_good_markers(inclidx)=1;
ts_vect=nan(length(trials_with_good_markers),1);
te_vect=nan(length(trials_with_good_markers),1);
ts_vect(logical(trials_with_good_markers))=time_start;
te_vect(logical(trials_with_good_markers))=time_end;

trials_used = strcmp({surf_str.out},'Success')' & strcmp({surf_str.type},'Normal')' & trials_with_good_markers;

if align==1
    ta_vect=ts_vect;
elseif align==2
    ta_vect=te_vect;
elseif align==3;
    [ta_vect,trials_used]=make_align_vector_from_BRT_between_start_and_end(fid,ts_vect,te_vect,trials_used);
else
    error('Align mode not supported')
end


max_align_start = max(ta_vect(trials_used)-ts_vect(trials_used));
align_offset = round(max_align_start/tres)+1;
max_align_end = max(te_vect(trials_used)-ta_vect(trials_used));
max_bind_width = align_offset + round(max_align_end/tres);
time_axis=linspace(min(ts_vect(trials_used)-ta_vect(trials_used)),max(te_vect(trials_used)-ta_vect(trials_used)),max_bind_width);
time_axis=time_axis(1:psthdec:end);


%keyboard
%%durations=te_vect-ts_vect;
%keyboard
%Use max(durations) to predefine the size of the nan matrix where psth will
%occur.
%%[~,trialID]=max(durations.*trials_used); %the max durations of the trials used

%%longest_psth = calc_psth(spks,ts_vect(trialID),te_vect(trialID),ts_vect(trialID),tres,psthfilter.Num,1,psthfilter.delay,psthdec);

%%grand_psth_matrix=nan(length(trials_used),length(longest_psth.data));
grand_psth_matrix=nan(length(trials_used),length(time_axis));


max_bins_to_delete=floor(psthfilter.length_ms/1e3/(tres/1e9*psthdec));
for trialID=1:length(trials_used)

    if trials_used(trialID)
        psth=calc_psth(spks,ts_vect(trialID),te_vect(trialID),ts_vect(trialID),tres,psthfilter.Num,1,psthfilter.delay,psthdec);
        if length(psth.data)>=max_bins_to_delete
            psth.data(1:max_bins_to_delete)=nan;
            psth.data(end-(max_bins_to_delete-1):end)=nan;
        else
            psth.data=nan(size(psth.data));
        end
        temp=size(grand_psth_matrix,2);
        if align==1
            % t_cac=psthfilter.% remove the 100ms at the beginning and the end of the psth
            grand_psth_matrix(trialID,1:length(psth.data))=psth.data;
        elseif align==2
            grand_psth_matrix(trialID,end-(length(psth.data)-1):end)=psth.data;
        elseif align==3
            [~,bin_start]=min(abs((time_axis-(ts_vect(trialID)-ta_vect(trialID)))));
            bind_end=bin_start+length(psth.data)-1;
            grand_psth_matrix(trialID,bin_start:bind_end)=psth.data;
        else
            error('Alignment mode not supported')
        end
        if size(grand_psth_matrix,2)~=temp
            time_axis(end+1)=time_axis(end)+time_axis(2)-time_axis(1);
        end
        
    end
end

time_axis=time_axis/1e9; %converted from ns to seconds
if align==1
    time_axis=time_axis-unique(surf_str_extra.pausedur(trials_used))+start_trig.offset;
elseif align==2
    time_axis=time_axis+end_trig.offset;
end


function [ta_vect,trials_used]=make_align_vector_from_BRT_between_start_and_end(fid,ts_vect,te_vect,trials_used)

ta_vect=nan(size(ts_vect));
ta_vect(trials_used)=double(0);
brt_onset_times = dhfun('GETMARKER',fid,'BRT_target_onset');

for i=1:length(ts_vect)
    if ~isnan(ta_vect(i))
        temp=brt_onset_times(brt_onset_times>=ts_vect(i) & brt_onset_times<=te_vect(i));
        if length(temp)>1
            error(['No exact match found for BRT inside trial number ' num2str(i)]);
        elseif length(temp)==1
            ta_vect(i)=temp;
        else
            ta_vect(i)=nan;
            trials_used(i)=false;
        end
    end
    
end

