function time_axis=calculate_time_axis_group_GLM(all_times,mat_used)

if mat_used<3
[~,Im]=max(cellfun(@length,all_times));
time_axis=all_times{Im};
elseif mat_used==3
    start_points= cellfun(@(x)x(1),all_times);
    end_points= cellfun(@(x)x(end),all_times);
    t_start=min(start_points);
    t_end=max(end_points);
    [~,Im]=max(cellfun(@length,all_times));
    tres=mode(diff(all_times{Im}));
    time_axis=t_start:tres:t_end;
end
    
    