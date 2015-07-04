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
    a=[total_behav.sig];
    
    tt(1) = all([a.discrimination_target]<(0.05/length(a)));
    tt(2) = all([a.detection_target]<(0.05/length(a)));
    tt(3) = all([a.discrimination_distractor]>(0.05/length(a)));
    tt(4) = all([a.detection_distractor]>(0.05/length(a)));
    
    disp(['Is all okay with monkey ' num2str(mm) '? ' num2str(all(tt))]) 
    mean_discr(mm,:) = [mean([b.discrimination_target]) mean([b.discrimination_distractor])];
    sem_discr(mm,:) = [std([b.discrimination_target]) std([b.discrimination_distractor])];%/sqrt(length(b)-1);
    
    mean_detect(mm,:) = [mean([b.detection_target]) mean([b.detection_distractor])];
    sem_detect(mm,:) = [std([b.detection_target]) std([b.detection_distractor])];%/sqrt(length(b)-1);


    mean_overall(mm,:) = [mean([b.detection_target].*[b.discrimination_target]) mean([b.detection_distractor].*[b.discrimination_distractor])]/100;
    sem_overall(mm,:) = [std([b.detection_target].*[b.discrimination_target]) std([b.detection_distractor].*[b.discrimination_distractor])]/100;%/sqrt(length(b)-1);
    
    overall_performance(mm,:) = mean_detect(mm,:) .* mean_discr(mm,:) * 1e-2;
    
    c=[total_behav.chance];
    chance_levels_detect_median(mm,:) = [mean([c.detection_target_median]) mean([c.detection_dist_median])];
    d1=[c.detection_target_confidence]; d2=[c.detection_dist_confidence];
    is_even=~logical(mod(1:length(d1),2));
    chance_levels_detect_confidence(mm,:) = [mean(d1(~is_even)) mean(d1(is_even)) mean(d2(~is_even)) mean(d2(is_even))];
    d=[c.discrimination_confidence];
    %chance_levels_discim_confidence(mm,:) = [50-mean(50-d(~is_even))/sqrt(length(b)-1) mean(d(is_even)-50)/sqrt(length(b)-1)+50];
    chance_levels_discim_confidence(mm,:) = [50-mean(50-d(~is_even)) mean(d(is_even)-50)+50];
    clear total_behav
end

f=figure;
set(f,'position',[440   328   1000   420])

subplot(1,3,1)
bh{1}=bar(mean_detect); hold on
ah(1)=gca;

subplot(1,3,2)
bh{2}=bar(mean_discr); hold on
ah(2)=gca;

subplot(1,3,3)
bh{3}=bar(mean_overall); hold on
ah(3)=gca;

titles = {'Detection','Discrimination', 'Overall'};
colors_used = distinguishable_colors(2);



for i=1:length(ah)
    set(ah(i),'ylim',[0 100],'Ytick',0:25:100,'Box','off','XtickLabel',{'M1','M2'}...
        );
    set(get(ah(i),'Ylabel'),'String','Performance (%)')
    set(get(ah(i),'Title'),'String',titles{i})
    set(bh{i}(1),'FaceColor',colors_used(1,:))
    set(bh{i}(2),'FaceColor',colors_used(2,:))
    bp=get(ah(i),'children');
    for j=1:length(bp)
        if i==1; used_error=sem_detect; elseif i==2; used_error=sem_discr; elseif i==3; used_error=sem_overall; else error('shit'); end
        x=get(get(bp(length(bp)+1-j),'children'),'xData');
        x=x(1,:)+(x(3,:)-x(2,:))/2;
        y=get(bp(length(bp)+1-j),'YData');
        subplot(1,length(ah),i)
%        errorbar(x,y,used_error(:,j),used_error(:,j),'color','k','linestyle','none','LineWidth',1)
    end
    if i==1
        for j=1:length(bp) %j=1:target; j=2:distractor
            x=get(get(bp(length(bp)+1-j),'children'),'xData');
            x=unique(x); % start and end coords of bars for M1 and M2
            Dx=(x(2)-x(1))/10;
            for mm=1:2
            x_coords=[x(1+2*(mm-1))-Dx x(2+2*(mm-1))+Dx];
            line(x_coords,[chance_levels_detect_median(mm,j) chance_levels_detect_median(mm,j)],'color','k')
            line(x_coords,[chance_levels_detect_confidence(mm,1+(j-1)*2) chance_levels_detect_confidence(mm,1+(j-1)*2)],'color','k','linestyle','--')
            line(x_coords,[chance_levels_detect_confidence(mm,2+(j-1)*2) chance_levels_detect_confidence(mm,2+(j-1)*2)],'color','k','linestyle','--')
            end
        end
    elseif i==2
        lh=line(xlim,[50 50],'color','k');
        tt=mean(chance_levels_discim_confidence);
        for j=1:2
            lh=line(xlim,[tt(j) tt(j)],'color','k','linestyle','--');
        end
    elseif i==3
        lgh = legend(ah(i),{'Cued Surface','Distractor Surface'});
        set(lgh,'Box','off')
    end
end

if save_fig
addpath('../export_fig/')
export_fig /Freiwald/ppolosecki/harbor/performance_figure.eps -eps -transparent
end