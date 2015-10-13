function [f]=plot_GLM_contrasts_for_presentation(results,true_contrasts,use_BRTs,cell_str,cell_no)

%INPUTS:
%   results: results structure as provided by make_GLM_and_contrasts_from_inst_firing

is_zscored=false;
no_error_bars=true;
if use_BRTs
    num_matrices=3;
    subplot_column_used=[1 3 2];
else
    num_matrices=2;
    subplot_column_used=[1 2];
end

show_legends=0;
lims_used=zeros(size(true_contrasts,1),num_matrices,2);
y_cell=cell(1,num_matrices);
y_std_cell=cell(1,num_matrices);

colors_used=distinguishable_colors(size(true_contrasts,1));

if size(true_contrasts,1)>4
  colors_used(4:end,:)=circshift(colors_used(4:end,:),[1 0]);
end



%colors_called=[1 2 4 7 5 3 6];
colors_called=[1 2 4 7 5 3 8 6];
if length(colors_called)<size(true_contrasts,1)
    error('Please, included more colors in variable colors_called')
end

% if size(true_contrasts,1)>=6
%     colors_used(1:6,:)=colors_used([1 2 6 4 5 3],:);
% end

contrasts_plotted={logical([0 0 0 0 0 0 0 0 0]);
    logical([0 0 0 0 0 0 0 0 0]);
    logical([0 0 0 0 0 0])};


for i=1:size(true_contrasts,1)
    f(i)=figure;
    %    for kk=size(true_contrasts,1):-1:i
    %    set(0,'currentfigure',f(i))
    
    set(f(i),'Position',get(0,'ScreenSize'));
    
    
    for mat_used=1:num_matrices
        subplot(1,num_matrices,subplot_column_used(mat_used)); hold on
        if mat_used==1
            xlim([-0.5 1.5]);%xlim([-1.5 1.5])
        elseif mat_used==2
            xlim([-0.8 0.025]);%0.5])
        elseif mat_used==3
            xlim([-.4 .4])
        else
            error('what mat??')
        end
        
        contrasts_plotted{true_contrasts(i,1)}(true_contrasts(i,2))=1; %note how this updates constrasts one by one
        
        y_matrix=[results{mat_used}.GLM(1).ces(contrasts_plotted{1},:);results{mat_used}.GLM(2).ces(contrasts_plotted{2},:);results{mat_used}.GLM(3).ces(contrasts_plotted{3},:)];
        y__std_matrix=[results{mat_used}.GLM(1).ces_std(contrasts_plotted{1},:);results{mat_used}.GLM(2).ces_std(contrasts_plotted{2},:);results{mat_used}.GLM(3).ces_std(contrasts_plotted{3},:)];
        if no_error_bars
            set(gca,'ColorOrder',colors_used(sort(colors_called(1:i)),:))
            temp=plot(results{mat_used}.time',y_matrix');
            set(temp,'linewidth',5)
            
        else
            [temp]= shadowcaster_ver3PP(results{mat_used}.time',y_matrix', y__std_matrix', [],colors_used(sort(colors_called(1:i)),:));
        end
    end
    if mat_used==1
        linehandles=temp;
    end
    clear temp;
    lims_used(i,mat_used,:)=ylim;
    %    end
end
mins=lims_used(:,:,1);
maxs=lims_used(:,:,2);
for i=1:size(true_contrasts,1)
    set(0,'currentfigure',f(i));
    for mat_used=1:num_matrices
        subplot(1,num_matrices,subplot_column_used(mat_used))
        pbaspect('manual')
        pbaspect([207 107 1])
        set(gca,'FontSize',25)
        
        ylim([min(mins(:)) max(maxs(:))])
        line(xlim,[0 0],'Color','black','LineStyle','-')
        line([0 0],ylim,'Color','black','LineStyle','--')
        xlabel ('Time from trigger event(s)')
        if mat_used==1
            xl=xlim;
            line([xl(1) xl(1)]+0.005,ylim,'Color','black','LineStyle','-','LineWidth',.5)
            if is_zscored
            ylabel('Effect Size (z-score units)')
            else
            ylabel('Effect Size (Hz)')
            end
            line([-1 -1],ylim,'Color',[.3 .3 .3],'LineStyle',':');
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
end

if nargin>3
    [axh,labelh]=suplabel([cell_str(cell_no).monkey '__' cell_str(cell_no).area '__cell_no:__' sprintf('%02.0f',cell_no)],'t');
    set(labelh,'Interpreter','none');
end

