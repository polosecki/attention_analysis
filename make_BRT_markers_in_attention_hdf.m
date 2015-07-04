function make_BRT_markers_in_attention_hdf(synched_hdf_file,username)

%Open file
fid = dhfun('open',synched_hdf_file,'r+');
[OPNAMES,OPINFOS] = dhfun('GETOPERATIONINFOS',fid);

%If it has BRT markers, close it
if any(cellfun(@isfield,OPINFOS,repmat({'marker_BRT'},size(OPINFOS))))
    dhfun('close',fid)
else
    %Make a backup copy of the file
    [pathstr,filename,ext] = fileparts(synched_hdf_file);    
    copyfile(synched_hdf_file,fullfile(pathstr,[filename '_not_brt' ext]));
    
    %Define event bits
    event.bit_target = 1;
    event.brt_target = 2;
    event.bit_distractor = 3;
    event.brt_distractor = 4;
    event.vsync = 5;
    event.fixation = 6;
    event.trialstart = 7;
    event.juice = 8;
    
    
    n_events = dhfun('GETEV2SIZE',fid);
    [evts,evev] = dhfun('READEV2',fid,1,n_events);
    
    idx_brt = (evev==event.brt_target);
    idx_brt_offset = (evev==-event.brt_target);
    idx_brt_distractor = (evev==event.brt_distractor);
    idx_brt_distractor_offset = (evev==-event.brt_distractor);
    
    loop_indexes=zeros(size(evts));
    
    loop_indexes(idx_brt)=1;
    loop_indexes(idx_brt_offset)=2;
    loop_indexes(idx_brt_distractor)=3;
    loop_indexes(idx_brt_distractor_offset)=4;
    
    brt_marker_names={'BRT_target_onset','BRT_target_offset','BRT_distractor_onset','BRT_distractor_offset'};
    
    for named_indx=1:length(brt_marker_names)
        dhfun('SETMARKER',fid,brt_marker_names{named_indx},evts(loop_indexes==named_indx));
    end
    
    
    
    operation_info.marker_BRT='Yes';
    dhfun('close',fid,'Creation of BRT markers',username,'make_BRT_markers_in_attention_hdf',operation_info);
end
