function [f]=plot_GLM_contrasts_for_presentation(results,true_contrasts)
%INPUTS:
%   results: results structure as provided by make_GLM_and_contrasts_from_inst_firing
%   contrasts_plotted: 2-by-1 cell with logical vectors indicating the contrasts to be plotted from

show_legends=0;

lims_used=zeros(size(true_contrasts,1),2,2);
%colors_used=distinguishable_colors(sum(contrasts_plotted{1}+sum(contrasts_plotted{2})));
colors_used=distinguishable_colors(6);
colors_used=colors_used([1 2 6 4 5 3],:);
y_cell=cell(1,2);
y_std_cell=cell(1,2);
for i=1:size(true_contrasts,1)
    f(i)=figure;
    set(f(i),'Position',get(0,'ScreenSize'));
    for mat_used=1:2
        subplot(1,2,mat_used); hold on
        if mat_used==1
            xlim([-1.5 2])
        elseif mat_used==2
            xlim([-2 0.5])
        else
            error('what mat??')
        end
        
        contrasts_plotted={logical([0 0 0 0 0 0 0 0 0]);
            logical([0 0 0 0 0 0 0 0 0])};
        contrasts_plotted{true_contrasts(i,1)}(true_contrasts(i,2))=1;
        y_cell{mat_used} =[y_cell{mat_used}; results{mat_used}.GLM1.ces(contrasts_plotted{1},:);results{mat_used}.GLM2.ces(contrasts_plotted{2},:)];
        y_std_cell{mat_used}=[y_std_cell{mat_used}; results{mat_used}.GLM1.ces_std(contrasts_plotted{1},:);results{mat_used}.GLM2.ces_std(contrasts_plotted{2},:)];
        [temp]= shadowcaster_ver3PP(results{mat_used}.time',y_cell{mat_used}', 2.*y_std_cell{mat_used}', [],colors_used(1:size(y_cell{mat_used},1),:));
        if mat_used==1
            linehandles=temp;
        end
        clear temp;
        lims_used(i,mat_used,:)=ylim;
    end
end

mins=lims_used(:,:,1);
maxs=lims_used(:,:,2);
for i=1:size(true_contrasts,1)
    set(0,'currentfigure',f(i));
    for mat_used=1:2
        
        subplot(1,2,mat_used)
        ylim([min(mins(:)) max(maxs(:))]);
        line(xlim,[0 0],'Color','black','LineStyle','-')
        line([0 0],ylim,'Color','black','LineStyle','--')
        xlabel ('Time from trigger event(s)')
    if mat_used==1                
        xl=xlim;
        line([xl(1) xl(1)]+0.005,ylim,'Color','black','LineStyle','-','LineWidth',.5)
        ylabel('Effect Size (Hz)')
        line([-1 -1],ylim,'Color',[.3 .3 .3],'LineStyle',':');
        if show_legends
            lg=legend(linehandles,{results{1}.GLM1.contrast.name{contrasts_plotted{1}} results{1}.GLM2.contrast.name{contrasts_plotted{2}}});
            set(lg,'Location','SouthWest');
            set(lg,'Box','off')
            %         set(lg,'units','pixels');
            %         lp=get(lg,'outerposition');
            %         set(lg,'outerposition',[lp(1:2),50,lp(4)]);
        end
    end
end
end