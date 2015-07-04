function out_RTs = calculate_RTs_from_hdf(fid,tmap)

start_trig.type ='Marker';
start_trig.name ='Stimulus_onset';
start_trig.offset = 0;
end_trig.type ='Marker';
end_trig.offset =0.1;
end_trig.name ='Fixation_out';
end_trig.offset =-.5;



[time_start,time_end,inclidx] = find_ref_points_heiko_PP(fid,tmap,start_trig,end_trig);


