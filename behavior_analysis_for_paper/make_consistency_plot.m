close all; clear all
monkeys = {'Quincy','Michel'};
areas = {'PITd','LIP'};
save_fig=true;
addpath('..')
for mm=1:2
    for aa=1:2
        load([monkeys{mm} '_' areas{aa} '_performance.mat'])
        t=[repmat({monkeys{mm}},size(behavior_performance))];
        [behavior_performance.monkey] = t{:};
        t=[repmat({areas{aa}},size(behavior_performance))];
        [behavior_performance.area] = t{:};
        
        if ~exist('total_behav','var')
            total_behav = behavior_performance;
        else
            total_behav = [total_behav behavior_performance];
        end
        clear behavior_performance
    end
    b=[total_behav.performance];
    k=[total_behav.consistency];
    p_val=[k.p_val]; succcess=[k.success]; chance=[k.chance]; confid=[k.confidence];
    
    %Holm-Bonferroni corrected inference
    tt(1) = all(frmrHolmBonferoni([p_val.targ])<.05);
    tt(2) = all(frmrHolmBonferoni([p_val.both])<.05);
    tt(3) = all(frmrHolmBonferoni([p_val.dist])<.05);
    tt(4) = all(frmrHolmBonferoni([p_val.none])<.05);
    
    disp(['Is all okay with monkey ' num2str(mm) '? ' num2str(tt)]) 
    mean_success(mm,:) = [mean([succcess.targ]) mean([succcess.both]) mean([succcess.dist]) mean([succcess.none])];
    sem_success(mm,:) = [std([succcess.targ]) std([succcess.both]) std([succcess.dist]) std([succcess.none])];

    chance_levels_mean(mm,:) = [mean([chance.targ]) mean([chance.both]) mean([chance.dist]) mean([chance.none])];
    a1=[confid.targ]; a2=[confid.both];a3=[confid.dist];a4=[confid.none];
    is_even=~logical(mod(1:length(a1),2));
    chance_levels_confidence(mm,:) = [mean(a1(~is_even)) mean(a1(is_even)) ...
        mean(a2(~is_even)) mean(a2(is_even)) mean(a3(~is_even)) mean(a3(is_even)) ...
        mean(a4(~is_even)) mean(a4(is_even))];

    
   clear total_behav
end
%%
f=figure;
%set(f,'position',[440   328   1000   420]);
set(f,'Position',get(0,'ScreenSize'));
bh=bar(mean_success); hold on;
ah=gca;


colors_used = distinguishable_colors(4);
w=.3;
colors_used = w*colors_used+(1-w)*ones(size(colors_used));
colors_used(2:3,:)=colors_used(3:-1:2,:);


    set(ah,'ylim',[0 75],'Ytick',0:25:75,'Box','off','XtickLabel',{'M1','M2'}...
        );
    set(get(ah,'Ylabel'),'String','% trials ')
    %set(get(ah(i),'Title'),'String',titles{i})
    for b=1:length(bh)
    set(bh(b),'FaceColor',colors_used(b,:));
    set(bh(b),'EdgeColor',colors_used(b,:));
    end
    
    bp=get(ah,'children');
    for j=1:length(bp)
        x=get(get(bp(length(bp)+1-j),'children'),'xData');
        x=x(1,:)+(x(3,:)-x(2,:))/2;
        y=get(bp(length(bp)+1-j),'YData');
       errorbar(x,y,sem_success(:,j),sem_success(:,j),'color','k','linestyle','none','LineWidth',1)
    end
    for j=1:length(bp) %j=1:target; j=2:distractor
            x=get(get(bp(length(bp)+1-j),'children'),'xData');
            x=unique(x); % start and end coords of bars for M1 and M2
            Dx=(x(2)-x(1))/10;
            for mm=1:2
            x_coords=[x(1+2*(mm-1))-Dx x(2+2*(mm-1))+Dx];
            line(x_coords,[chance_levels_mean(mm,j) chance_levels_mean(mm,j)],'color','k')
            line(x_coords,[chance_levels_confidence(mm,1+(j-1)*2) chance_levels_confidence(mm,1+(j-1)*2)],'color','k','linestyle','--')
            line(x_coords,[chance_levels_confidence(mm,2+(j-1)*2) chance_levels_confidence(mm,2+(j-1)*2)],'color','k','linestyle','--')
            end
        end
        

        lgh = legend(ah,{'Target','Both','Distractor','None'});
        set(lgh,'Box','off')
        title('Trial outcomes sorted by consitency with cued and distractor MDS')

if save_fig
addpath('../export_fig/')
export_fig /Freiwald/ppolosecki/harbor/consistency_figure.eps -eps -transparent
end