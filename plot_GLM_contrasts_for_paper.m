function [linehandles]=plot_GLM_contrasts_for_paper(results,contrasts_plotted,use_BRTs,cell_str,cell_no)
%INPUTS:
%   results: results structure as provided by make_GLM_and_contrasts_from_inst_firing
%   contrasts_plotted: 2-by-1 cell with logical vectors indicating the contrasts to be plotted from

is_zscored=true;%true;
use_shadowcaster=false;%true;
if use_BRTs
    num_matrices=3;
    subplot_column_used=[1 3 2];
else
    num_matrices=2;
    subplot_column_used=[1 2];
end

show_legends=0;
f=figure;
set(f,'Position',get(0,'ScreenSize'));
lims_used=zeros(num_matrices,2);
%colors_used=distinguishable_colors(sum(contrasts_plotted{1})+sum(contrasts_plotted{2})+sum(contrasts_plotted{3}));
colors_used=distinguishable_colors(8);
%whos colors_used
colors_used(4:end,:)=circshift(colors_used(4:end,:),[1 0]);
%colors_used=distinguishable_colors(6);
%colors_used=colors_used([1 2 6 4 3 5],:);
if all(find(cellfun(@sum,contrasts_plotted))==3)
    colors_used=colors_used([1 2 7 5 4 8 3 6],:);
elseif any(find(cellfun(@sum,contrasts_plotted))==1)
    colors_used=colors_used([7 5 1 2 4 8 3 6],:);
end

%% DELETE THIS IF PLOTTING MRE THAN 3 LINES:
if sum(cellfun(@sum,contrasts_plotted))<=3
colors_used=colors_used([3,4,8],:);
end

%%
%colors_used=distinguishable_colors(6);
%colors_used=colors_used([1 2 6 4 3 5],:);
for mat_used=1:num_matrices
    subplot(1,num_matrices,subplot_column_used(mat_used)); hold on

    if mat_used==1
        %xlim([-1.5 1.5])
        xlim([-.5 1.5])
    elseif mat_used==2
%        xlim([-0.6 0.025]);%0.5])
        xlim([-0.8 0.025]);%0.5])

    elseif mat_used==3
        xlim([-.4 .4])
    else
        error('what mat??')
    end

    
    y_matrix=[results{mat_used}.GLM(1).ces(contrasts_plotted{1},:);results{mat_used}.GLM(2).ces(contrasts_plotted{2},:);results{mat_used}.GLM(3).ces(contrasts_plotted{3},:)];
    y__std_matrix=[results{mat_used}.GLM(1).ces_std(contrasts_plotted{1},:);results{mat_used}.GLM(2).ces_std(contrasts_plotted{2},:);results{mat_used}.GLM(3).ces_std(contrasts_plotted{3},:)];
    if use_shadowcaster        
    [temp]= shadowcaster_ver3PP(results{mat_used}.time',y_matrix', y__std_matrix', [],colors_used);
        set(temp,'linewidth',5)

    else
        set(gca,'ColorOrder',colors_used)
        temp=plot(results{mat_used}.time',y_matrix');
        set(temp,'linewidth',5)
        
    end
    %2*y__std_matrix'
    if mat_used==1
        linehandles=temp;
    end
    clear temp;
    lims_used(mat_used,:)=ylim;
end

for mat_used=1:num_matrices
    subplot(1,num_matrices,subplot_column_used(mat_used))
    pbaspect('manual')
    pbaspect([1 1 1])
%    pbaspect([207 107 1])
    set(gca,'FontSize',25)
    if ~is_zscored
        ylim([min(lims_used(:,1)) max(lims_used(:,2))]);
    else
    ylim([-1.2 2.2])
%    ylim([-1.2 3])
%    ylim([-.6 2.2])
    %ylim([-40 60])
    yt=get(gca,'ytick');
    set(gca,'ytick',yt(rem(yt,1)==0));
        end
    line(xlim,[0 0],'Color','black','LineStyle','-')
    line([0 0],ylim,'Color','black','LineStyle','--')
    if mat_used==1
        xlabel ('Time from stimulus onset (s)')
    elseif mat_used==2
        xlabel ('Time from saccade onset (s)')
    elseif mat_used==3
        xlabel ('Time from BRT onset (s)')
    end
    ylabel('Effect Size (z-score)')
    
    if mat_used==1
        
        if min(xlim)<-1
        line([-1 -1],ylim,'Color',[.3 .3 .3],'LineStyle',':');
        end
        if show_legends
%             a={results{1}.GLM(1).contrast.name{contrasts_plotted{1}} results{1}.GLM(2).contrast.name{contrasts_plotted{2}}};
%             for hh=1:length(a)
%                 disp([num2str(hh) ' ' a{hh}]);
%             end
        lg=legend(linehandles,{results{1}.GLM(1).contrast.name{contrasts_plotted{1}} results{1}.GLM(2).contrast.name{contrasts_plotted{2}} results{1}.GLM(3).contrast.name{contrasts_plotted{3}}});
        set(lg,'Location','SouthWest');
        set(lg,'Box','off')
        %         set(lg,'units','pixels');
        %         lp=get(lg,'outerposition');
        %         set(lg,'outerposition',[lp(1:2),50,lp(4)]);
        end
    end

end


if nargin>3 & ~isempty(cell_str)
[axh,labelh]=suplabel([cell_str(cell_no).monkey '__' cell_str(cell_no).area '__cell_no:__' sprintf('%02.0f',cell_no)],'t');
set(labelh,'Interpreter','none');
end

