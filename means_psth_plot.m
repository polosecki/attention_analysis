function []=means_psth_plot(mean_mat,std_mat,plot_params)

hold all
t=plot_params('t');


[linehandle]= shadowcaster_ver3PP(t', mean_mat', std_mat', [],plot_params('colors'));

if ismember('linestyle',plot_params.keys) && ~isempty(plot_params('linestyle'))
    ls=plot_params('linestyle');
    for i=1:length(linehandle)
    set(linehandle(i),'LineStyle',ls{i});
    end
end

xlim(plot_params('xlims'))
ylim(plot_params('ylims'))
pbaspect(plot_params('pbaspect'))
box off
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--')    
