clc;clear all; close all
monkey={'Quincy','Michel'};
area={'PITd','LIP'};

use_high_res_data=1;
redo_stats=0;
use_mean_firing=0;
remove_mgs_saccade_epoch=1; %This might be desirable, since we already have its difference with mem epoch as a variable 
save_figures_to_harbor=0;
figure_dir='/Freiwald/ppolosecki/harbor';
noise_model='poisson';
%% Make or load statistics for each area and monkey
addpath('../')
if redo_stats
    for mm=1:length(monkey)
        for aa=1:length(area)
            [cell_ces, cell_sig, param_names]=calculate_stats_from_cell_list(monkey{mm},area{aa},noise_model);
            [timing]= measure_timing_of_effects_fun(monkey{mm},area{aa},use_high_res_data);
            [RT_data]= get_RT_from_cell_list(monkey{mm},area{aa});
            full_stats(aa,mm).ces=cell_ces;
            full_stats(aa,mm).sig=cell_sig;
            full_stats(aa,mm).param_names=param_names;
            full_stats(aa,mm).monkey=monkey{mm};
            full_stats(aa,mm).area=area{aa};
            full_stats(aa,mm).timing=timing;
            full_stats(aa,mm).RT_data=RT_data;
            save full_stats full_stats
        end
    end
    sendmail('9178422510@txt.att.net', 'Function analyses report', 'Loop completed');
else
    load full_stats
end

if ~use_mean_firing & any(strcmp(full_stats(1,1).param_names,'mean firing'))
    mf_index=find(strcmp(full_stats(1,1).param_names,'mean firing'));
    for mm=1:length(monkey)
        for aa=1:length(area)
            full_stats(aa,mm).ces(:,[mf_index mf_index+1])=[];
            full_stats(aa,mm).sig(:,[mf_index mf_index+1])=[];
            full_stats(aa,mm).param_names([mf_index mf_index+1])=[];
        end
    end
    clear mf_index
end

if remove_mgs_saccade_epoch
    mgs_sacc_index=find(strcmp(full_stats(1,1).param_names,'mgs saccade-epoch'));
    for mm=1:length(monkey)
        for aa=1:length(area)
            full_stats(aa,mm).ces(:,mgs_sacc_index)=[];
            full_stats(aa,mm).sig(:,mgs_sacc_index)=[];
            full_stats(aa,mm).param_names(mgs_sacc_index)=[];
        end
    end
    clear mgs_sacc_index
end

addpath('..')

%% Make plots of saccade lattencies

x_axis=-.8:.05:0;

f=figure; 
set(f,'Position',get(0,'ScreenSize')); 

for mm=1:2
    subplot(1,2,mm);
    hold all
    for aa=1:2
t_vect=full_stats(aa,mm).timing.saccade.onset;
[n]=hist(t_vect,x_axis);
plot(x_axis,cumsum(n)/length(t_vect)*100,'.-')
xlabel('Time before saccade (s)')
ylabel('Population fraction (%)')
title(monkey{mm})
end
lg=legend({full_stats(:,mm).area});
set(lg,'box','off')
end

[axh,labelh]=suplabel('Cummulative saccade latencies','t');
set(labelh,'Interpreter','none');
set(labelh,'FontSize',get(labelh,'FontSize')+1)
%% Make plot of attention latencies
x_axis=-1:.05:1;

f=figure; 
set(f,'Position',get(0,'ScreenSize')); 

for mm=1:2
    subplot(1,2,mm);
    hold all
    for aa=1:2
        cells_used=(abs(full_stats(aa,mm).timing.attention.duration)>0.05);
        t_vect=full_stats(aa,mm).timing.attention.onset;
        [n]=hist(t_vect(cells_used),x_axis);
        plot(x_axis,cumsum(n)/length(t_vect)*100,'.-')
        xlabel('Time from surface onset (s)')
        ylabel('Population fraction (%)')
        title(monkey{mm})
    end
lg=legend({full_stats(:,mm).area});
set(lg,'box','off')
end

[axh,labelh]=suplabel('Cummulative attention latencies','t');
set(labelh,'Interpreter','none');
set(labelh,'FontSize',get(labelh,'FontSize')+1)

%% Make statistical tests of differences in saccade lattencies
saccade_times=false(0);
sacc_is_pitd=false(0);
sacc_is_quincy=false(0);
for mm=1:2
    for aa=1:2
       temp=full_stats(aa,mm).timing.saccade.onset;
       good_idx=~isnan(temp);
       saccade_times=[saccade_times; temp(good_idx)];
       if strcmp(full_stats(aa,mm).area,'PITd')
           sacc_is_pitd=[sacc_is_pitd; true(sum(good_idx),1)];
       else
           sacc_is_pitd=[sacc_is_pitd; false(sum(good_idx),1)];
       end
       if strcmp(full_stats(aa,mm).monkey,'Quincy')
           sacc_is_quincy=[sacc_is_quincy; true(sum(good_idx),1)];
       else
           sacc_is_quincy=[sacc_is_quincy; false(sum(good_idx),1)];
       end
    end
    %The difference in saccade latencies seems to be not-significant:
    [p,h,stats]=ranksum(saccade_times(sacc_is_pitd),saccade_times(~sacc_is_pitd))
end
%both monkeys:
    [p,h,stats]=ranksum(saccade_times(sacc_is_pitd),saccade_times(~sacc_is_pitd))
%Quincy:
   [p,h,stats]=ranksum(saccade_times(sacc_is_pitd & sacc_is_quincy),saccade_times(~sacc_is_pitd & sacc_is_quincy))
%Michel:
   [p,h,stats]=ranksum(saccade_times(sacc_is_pitd & ~sacc_is_quincy),saccade_times(~sacc_is_pitd & ~sacc_is_quincy))

   
%The difference in saccade latencies seems to be not-significant:
%[p,h,stats]=ranksum(saccade_times(sacc_is_pitd),saccade_times(~sacc_is_pitd));%,'method','exact')



%% Makes a matrix por PCA input and vectors identifying areas



% [coeff,score,latent,tsquared,explained,mu] = pca(X,Name,Value)
%coeff: says the pca coefficients
%score: each column of score represents one principal component
%latent: the variance explained by each component
%Xcentered = score*coeff'
%use: 'Centered',false to avoid the mean of each column (i.e., neuron) from being removed
%Num components: 'NumComponents',12 % to get 12 PCs only

pca_input_matrix=[];
is_pitd=[];
is_monkey1=[];
for aa=1:2
    for mm=1:2
        idx=~isnan(full_stats(aa,mm).ces(:,1));
        if aa==1
            is_pitd=[is_pitd; true(sum(idx),1)];
        else
            is_pitd=[is_pitd; false(sum(idx),1)];
        end
        if mm==1
            is_monkey1=[is_monkey1; true(sum(idx),1)];
        else
            is_monkey1=[is_monkey1; false(sum(idx),1)];
        end
        pca_input_matrix=[pca_input_matrix; full_stats(aa,mm).ces(idx,:)];
    end
end

is_pitd=logical(is_pitd);
is_monkey1=logical(is_monkey1);
[coeff,score,latent,~,explained,mu] = pca(pca_input_matrix,'NumComponents',4);

pca_elements.pca_input_matrix=pca_input_matrix;
pca_elements.coeff=coeff;
pca_elements.score=score;
pca_elements.explained=explained;
pca_elements.mu=mu;
pca_elements.column_names=full_stats(1,1).param_names;

%% Plot inter-area distributions accross properties
addpath('../spaceplots/')

num_histnorm_cols=10;
f=figure;
set(f,'Position',get(0,'ScreenSize')); 
for par_num=1:size(pca_input_matrix,2)
    for mm=1:2
        subplot_index=(mm-1)*size(pca_input_matrix,2)+par_num;
        subplot(2,size(pca_input_matrix,2),subplot_index)
        if mm==1; monkey_used=is_monkey1;
        else monkey_used=~is_monkey1; end
        
        histnorm(pca_input_matrix(is_pitd & monkey_used,par_num),num_histnorm_cols)
        %histnorm(pca_input_matrix(is_pitd,4),15)
        set(gca,'color','none')
        h = findobj(gca,'type','patch');
        set(h,'FaceColor','r','EdgeColor','none','facealpha',0.7)
        hold on
        histnorm(pca_input_matrix(~is_pitd & monkey_used,par_num),num_histnorm_cols)
        %histnorm(pca_input_matrix(~is_pitd,4),15)
        h = findobj(gca,'Type','patch');
        set(h,'facealpha',0.7,'EdgeColor','none');
        ylabel('counts')
        xlabel('z-score')
        p=ranksum(pca_input_matrix(is_pitd & monkey_used,par_num),pca_input_matrix(~is_pitd & monkey_used,par_num));
        %[~,p] = ttest2(pca_input_matrix(is_pitd & monkey_used,par_num),pca_input_matrix(~is_pitd & monkey_used,par_num))
%         if subplot_index==1
%         legend('PITd','LIP','location','North');%'WestOutside')
%         end
        title(sprintf('%s p=%01.2g',full_stats(1,1).param_names{par_num},p),'fontsize',7);
    end
end
set(f,'position',[31         538        1405         234])
spaceplots(gcf,[0 0 0 0],[.02 .02]);
%tightfig
% 
num_histnorm_cols=15;
f=figure;
set(f,'Position',get(0,'ScreenSize')); 

for par_num=1:size(pca_input_matrix,2)
        subplot_index=par_num;
        subplot(1,size(pca_input_matrix,2),subplot_index)        
        histnorm(pca_input_matrix(is_pitd,par_num),num_histnorm_cols)
        %histnorm(pca_input_matrix(is_pitd,4),15)
        set(gca,'color','none')
        h = findobj(gca,'type','patch');
        set(h,'FaceColor','r','EdgeColor','none','facealpha',0.7)
        hold on
        histnorm(pca_input_matrix(~is_pitd,par_num),num_histnorm_cols)
        %histnorm(pca_input_matrix(~is_pitd,4),15)
        h = findobj(gca,'Type','patch');
        set(h,'facealpha',0.7,'EdgeColor','none');
        ylabel('counts')
        xlabel('z-score')
%         if subplot_index==1
%         legend('PITd','LIP','location','North')
%         end
 p=ranksum(pca_input_matrix(is_pitd,par_num),pca_input_matrix(~is_pitd,par_num));
 %[~,p]=ttest2(pca_input_matrix(is_pitd,par_num),pca_input_matrix(~is_pitd,par_num))
        title(sprintf('%s p=%01.2g',full_stats(1,1).param_names{par_num},p),'fontsize',7);
end
set(f,'position',[31         652        1405         120])
spaceplots(gcf,[0 0.05 0 0],[.02 .02])



%% Plot correlation between statistics for each monkey and area:

f=figure;
set(f,'Position',get(0,'ScreenSize')); 

load('../CoolWarmMap.mat')
cmap(1,:)=[];
cmap(end,:)=[];

for mm=1:length(monkey)
    for aa=1:length(area)
        subplot(2,2,length(monkey)*(aa-1)+mm)
        idx=~isnan(full_stats(aa,mm).ces(:,1));
        corr_mat{aa,mm}=corr(full_stats(aa,mm).ces(idx,:));
        imagesc(corr_mat{aa,mm},[-1 1])
        t=title([area{aa} '\_' monkey{mm}]);
        colormap(cmap)
        set(gca,'ytick',1:length(full_stats(1,1).param_names))
        set(gca,'yticklabel',full_stats(1,1).param_names)
        set(gca,'xtick',1:length(full_stats(1,1).param_names))
set(gca,'xticklabel',full_stats(1,1).param_names)
hh = rotateXLabels(gca,20,'MaxStringLength',25);
    end
end
hp=get(subplot(2,2,4),'Position')
colorbar('Position',[hp(1)+hp(3)+0.02 hp(2), 0.02 hp(2)+hp(3)*2.1])
%f0=myaa(2)
if save_figures_to_harbor
    filename=['population_stats_figure_' sprintf('%02.0f',f0)];
    saveas(f0,fullfile(figure_dir,filename),'png');
    print(f0,'-dpng',fullfile(figure_dir,[filename '.png']))
end

%%
% Plot the same, collapsed across monkeys
f=figure;
set(f,'Position',get(0,'ScreenSize'));
for aa=1:2
    if aa==1; i1=is_pitd;else; i1=~is_pitd;end
subplot(1,2,aa);
general_corr_matrix=corr(pca_input_matrix(i1,:));
imagesc(general_corr_matrix,[-1 1])
colormap(cmap)

set(gca,'ytick',1:length(full_stats(1,1).param_names))
        set(gca,'yticklabel',full_stats(1,1).param_names)
        set(gca,'xtick',1:length(full_stats(1,1).param_names))
  set(gca,'xticklabel',full_stats(1,1).param_names)
 hh = rotateXLabels(gca,20,'MaxStringLength',25);
     
title(area{aa})
end

hp=get(subplot(1,2,2),'Position')
colorbar('Position',[hp(1)+hp(3)+0.02 hp(2), 0.02 hp(2)+hp(3)*2.1])

f0=myaa(2)
%close(f)
%%
%Plot the difference of correlations between areas  
f=figure;
set(f,'Position',get(0,'ScreenSize'));
%cc1=(corr_mat{1,1}+corr_mat{1,2})/2-(corr_mat{2,1}+corr_mat{2,2})/2;
cc1=corr(pca_input_matrix(is_pitd,:))-corr(pca_input_matrix(~is_pitd,:));
%general_corr_matrix=corr(pca_input_matrix);
%imagesc(general_corr_matrix,[-1 1])
imagesc(cc1,[-1 1])
colormap(cmap)
title('Mean correlation differences (PITd-LIP)')
set(gca,'yticklabel',full_stats(1,1).param_names)
set(gca,'xticklabel',full_stats(1,1).param_names)
 hh = rotateXLabels(gca,25);
%spaceplots(gcf,[0.08 0.05 0.05 0.2],[.02 .02])
colorbar
f0=myaa(2)
if save_figures_to_harbor
                filename=['population_stats_figure_' sprintf('%02.0f',f0)];
 saveas(f0,fullfile(figure_dir,filename),'png');
end
%close(f)

% Plot correlations in both monkeys, collapsed accross areas
f=figure;
set(f,'Position',get(0,'ScreenSize')); 
%cc2=((corr_mat{1,1}+corr_mat{1,2})/2+(corr_mat{2,1}+corr_mat{2,2})/2)/2;
%imagesc(cc2,[-1 1])
general_corr_matrix=corr(pca_input_matrix);
imagesc(general_corr_matrix,[-1 1])

colormap(cmap)
title('Overall correlations (PITd andLIP)')
set(gca,'yticklabel',full_stats(1,1).param_names)
set(gca,'yticklabel',full_stats(1,1).param_names)
set(gca,'xtick',get(gca,'ytick'))
set(gca,'xticklabel',full_stats(1,1).param_names)
hh = rotateXLabels(gca,25);
 
%spaceplots(gcf,[0.08 0.05 0.05 0.2],[.02 .02])
colorbar
f0=myaa(2)
if save_figures_to_harbor
                filename=['population_stats_figure_' sprintf('%02.0f',f0)];
 saveas(f0,fullfile(figure_dir,filename),'png');
end
%close(f)

%% Plot variables of interes against each other for both monkeys and areas

x_ind=5;%4;
y_ind=7;%10;


max_m=2;
max_a=2;
f=figure;
set(f,'Position',get(0,'ScreenSize'));
for mm=1:max_m
    for aa=1:max_a
        subplot(max_a,max_m,max_a*(aa-1)+mm);
        plot(full_stats(aa,mm).ces(:,x_ind),full_stats(aa,mm).ces(:,y_ind),'.')
        axis equal
        lsline;
        ok_cells=~isnan(full_stats(aa,mm).ces(:,x_ind));
        [r,p]=corrcoef(full_stats(aa,mm).ces(ok_cells,x_ind),full_stats(aa,mm).ces(ok_cells,y_ind));
        title([full_stats(aa,mm).area ' ' full_stats(aa,mm).monkey ' r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
    end
end
        suplabel(full_stats(1,1).param_names{x_ind},'x')
        suplabel(full_stats(1,1).param_names{y_ind},'y');

% Plot the same, collapsed across monkeys
f=figure;
set(f,'Position',get(0,'ScreenSize'));
for aa=1:2
    if aa==1; i1=is_pitd;else; i1=~is_pitd;end
subplot(1,2,aa); hold all
plot(pca_input_matrix(i1,x_ind),pca_input_matrix(i1,y_ind),'.')
axis equal
lsline
[r,p]=corrcoef(pca_input_matrix(i1,x_ind),pca_input_matrix(i1,y_ind))
xlabel([pca_elements.column_names{x_ind} ' (z-score)'])
ylabel([pca_elements.column_names{y_ind} ' (z-score)']);
title([full_stats(aa,mm).area  ' r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
end
suplabel([pca_elements.column_names{y_ind} ' vs. ' pca_elements.column_names{x_ind}],'t')

%% See if PITd and LIP can be separated using SVM

%features_for_decoding=[2 3 4];
features_for_decoding=1:size(pca_input_matrix,2); %use all features
xdata=pca_input_matrix(:,features_for_decoding);
group=is_pitd;

max_perm=50;
%posta=0.3113;
perform_decoding=0;
calculate_null_distribution=0;

if perform_decoding
    %c = cvpartition(sum(is_monkey1),'leaveout');
    c = cvpartition(length(group),'leaveout');
    %%c = cvpartition(length(group),'kfold',round(.9*length(group)));
    posta=crossval('mcr',xdata,group,'Predfun',@svm_decoding,'partition',c);
    
    if calculate_null_distribution
        for i=1:max_perm
            disp(['running permutation number ' num2str(i)])
            c = cvpartition(length(group),'leaveout');
            mcr(i) = crossval('mcr',xdata,group(randperm(length(group))),'Predfun',@svm_decoding,'partition',c);
        end
        save crosval_general mcr
    else
        load crosval_general
    end
    p_value=mean(mcr<posta);

figure;
histnorm(1-mcr)
box off
line([1-posta 1-posta],ylim,'Color','black')
line([.5 .5],ylim,'Color','black','LineStyle','--')
lg=legend({'null distribution','real cell population'});
set(lg,'Box','off')
ylabel('Probability density')
xlabel('Decoding performance')
title('Decoding PITd vs LIP')
perform_decoding_per_monkey=0;
calculate_null_distribution_per_monkey=0;

if perform_decoding_per_monkey
    c = cvpartition(sum(~is_monkey1),'leaveout');
    %%c = cvpartition(length(group),'kfold',round(.9*length(group)));
    posta=crossval('mcr',xdata(~is_monkey1,:),group(~is_monkey1),'Predfun',@svm_decoding,'partition',c) 
    if calculate_null_distribution_per_monkey
        for i=1:max_perm
            disp(['running permutation number ' num2str(i)])
            c = cvpartition(sum(~is_monkey1),'leaveout');
            g_data=group(~is_monkey1);
            mcr(i) = crossval('mcr',xdata(~is_monkey1,:),g_data(randperm(length(g_data))),'Predfun',@svm_decoding,'partition',c);
        end
        save crosval_per_monkey mcr
    else
        load crosval_per_monkey
    end
    p_value=mean(mcr<posta);
end
end
%% PLOT PCA components
figure;
hold all
pcx=1;
pcy=2;
%pcz=3;
colors={'blue','red'};
legend_names={'PITd','LIP'}

%label_used=is_pitd;
label_used=is_monkey1;
for aa=1:2
    if aa==1;indx=label_used;
    else indx=~label_used; end
% plot(pca_input_matrix(indx,pcx),pca_input_matrix(indx,pcy),'.');
%ah=plot(score(indx,pcx),score(indx,pcy),'b.');
ah=plot(score(indx,pcx),score(indx,pcy),'b.');
%ah=plot3(score(indx,pcx),score(indx,pcy),score(indx,pcz),'.')
set(ah,'color',colors{aa})
end
legend(legend_names)
xlabel(['PC ' num2str(pcx)])
ylabel(['PC ' num2str(pcy)]);

%% Plots MGS distributions
sig_matrix=[];
for aa=1:2
    for mm=1:2
        idx=~isnan(full_stats(aa,mm).ces(:,1));
        %sig_matrix=[sig_matrix;abs(full_stats(aa,mm).sig(idx,:))<0.05];
        sig_matrix=[sig_matrix;full_stats(aa,mm).sig(idx,:)]; %changed;
    end
end

%f=figure;
xlim_used=[-3 4.5];
logical_vals=[true;false];
contrasts_plotted=[10;12];%[7;8];
num_histnorm_cols=15;
monkey_used=3;
mult_compare_method='fdr';%'fdr';'holm-bonferroni';'bonferroni';
p_threshold=0.05;

if monkey_used==1
    monkey_vector=is_monkey1;
elseif monkey_used==2;
    monkey_vector=~is_monkey1;
else
    monkey_vector=true(size(is_monkey1));
end
monkey_names={'Quincy','Michel','Both monkeys'};

bar_colors(1, :) = [255 0 0]/255; %(darker) red
bar_colors(2, :) = [0 0 128]/255; %navy/dark blue

for num_area=1:2
    f(num_area)=figure;
    
    for pp=1:length(contrasts_plotted)
        
        figure(f(num_area)); colormap(gray)
        
        %subplot(2,length(parnums),(num_row-1)*length(parnums)+pp)
        subplot(1,length(contrasts_plotted),pp);
        par_num=contrasts_plotted(pp);
        switch mult_compare_method
            case 'none'
                sig_vector=sig_matrix(:,par_num)<p_threshold;
            case 'fdr'
                relevant_idx=(is_pitd==logical_vals(num_area)) & monkey_vector;
                tested_pvals=sig_matrix(relevant_idx,par_num);
                %[h_tested, crit_p, ~]=fdr_bh(abs(tested_pvals),p_threshold);
                [h_tested, crit_p]=fdr_bky(abs(tested_pvals),p_threshold);
                sig_vector=false(size(relevant_idx));
                sig_vector(relevant_idx)=h_tested;
            case 'bonferroni'
                sig_vector=abs(sig_matrix(:,par_num))<(p_threshold/sum(is_pitd==logical_vals(num_area) & monkey_vector));
            case 'holm-bonferroni'
                relevant_idx=(is_pitd==logical_vals(num_area)) & monkey_vector;
                tested_pvals=abs(sig_matrix(relevant_idx,par_num));
                adjustedPVals = frmrHolmBonferoni(tested_pvals);
                sig_vector=false(size(relevant_idx));
                sig_vector(relevant_idx)=adjustedPVals<p_threshold;
            otherwise
                error('Invalid choice of multiple comparison control')
        end
        edges=linspace(min(pca_input_matrix(:,par_num)), max(pca_input_matrix(:,par_num)),num_histnorm_cols);
        counts=nan(length(logical_vals),length(edges));
        for i=1:length(logical_vals)
            n = histc(pca_input_matrix(sig_vector==logical_vals(i) & is_pitd==logical_vals(num_area) & monkey_vector,par_num),edges);
            counts(i,:)=n;
        end
        n_cells=sum(counts(:));        
        counts=counts/n_cells*100;
        h=bar(edges,counts','stacked');

        %colormap(cmap)
        this_cmap=colormap;
        this_cmap(1,:)=bar_colors(num_area, :);
        colormap(this_cmap); %freezeColors;
        v=version('-release');
        for i=1:length(h)
            if str2num(v(1:end-1))>=2015
                set(h(i),'EdgeColor',bar_colors(num_area, :))
            else
                set(get(h(i),'Children'),'EdgeColor',bar_colors(num_area, :))
                %set(get(h(i),'Children'),'FaceColor',bar_colors(num_row, :))
            end
        end
        xlim(xlim_used)
        pbaspect([172 142 1])
        box off
        title(sprintf('%01.2g%% of%3.0f cells are significant in %s',sum(counts(1,:)),n_cells,monkey_names{monkey_used}))
        if max(get(gca,'ytick'))==10
            set(gca,'ytick',[0 5 10])
        elseif max(get(gca,'ytick'))==20
            set(gca,'ytick',[0 10 20])
        end
        line([0 0],ylim,'Color','black','LineStyle','--');
        
    end
end

save_figures_to_harbor=true;
if save_figures_to_harbor
    for i=1:length(f)
    filename=['MGS_stats_' area{i} '_' monkey_names{monkey_used}];
    %plot2svg(fullfile(figure_dir,[filename '.svg']),f(i));
    saveas(f(i),fullfile(figure_dir,filename),'fig');
    saveas(f(i),fullfile(figure_dir,filename),'epsc');
    end
end

%% Make nice histograms in both monkeys
contrasts_plotted=[1 4 5 6 2 3 8 7];
%contrasts_plotted=1:8;
%contrasts_plotted=[1 4 5 2 3];

monkey_used=3;
contrast_names=full_stats(1,1).param_names;    
if monkey_used==1
    monkey_vector=is_monkey1;
elseif monkey_used==2;
    monkey_vector=~is_monkey1;
else
    monkey_vector=true(size(is_monkey1));
end
monkey_names={'Quincy','Michel','Both monkeys'};
bar_colors(1, :) = [255 0 0]/255; %(darker) red
bar_colors(2, :) = [0 0 128]/255; %navy/dark blue
mult_compare_method='fdr';%'holm-bonferroni';%'fdr';%'none';%'fdr';%'holm-bonferroni';%'bonferroni';%'bonferroni';'fdr';'none';
p_threshold=0.05;
sig_matrix=[];
for aa=1:2
    for mm=1:2
        idx=~isnan(full_stats(aa,mm).ces(:,1));
        sig_matrix=[sig_matrix;full_stats(aa,mm).sig(idx,:)]; %changed;
    end
end

logical_vals=[true;false];
num_histnorm_cols=15;
f=[];
for num_area=1:2
    f(num_area)=figure;
    
    for pp=1:length(contrasts_plotted)
        
        figure(f(num_area)); colormap(gray)
        
        %subplot(2,length(parnums),(num_row-1)*length(parnums)+pp)
        subplot(2,ceil(length(contrasts_plotted)/2),pp);
        par_num=contrasts_plotted(pp);
        switch mult_compare_method
            case 'none'
                sig_vector=sig_matrix(:,par_num)<p_threshold;
            case 'fdr'
                relevant_idx=(is_pitd==logical_vals(num_area)) & monkey_vector;
                tested_pvals=sig_matrix(relevant_idx,par_num);
                [h_tested, crit_p]=fdr_bky(abs(tested_pvals),p_threshold);
                %[h_tested, crit_p, ~]=fdr_bh(abs(tested_pvals),p_threshold);
                sig_vector=false(size(relevant_idx));
                sig_vector(relevant_idx)=h_tested;
            case 'bonferroni'
                sig_vector=abs(sig_matrix(:,par_num))<(p_threshold/sum(is_pitd==logical_vals(num_area) & monkey_vector));
            case 'holm-bonferroni'
                relevant_idx=(is_pitd==logical_vals(num_area)) & monkey_vector;
                tested_pvals=abs(sig_matrix(relevant_idx,par_num));
                adjustedPVals = frmrHolmBonferoni(tested_pvals);
                sig_vector=false(size(relevant_idx));
                sig_vector(relevant_idx)=adjustedPVals<p_threshold;
            otherwise
                error('Invalid choice of multiple comparison control')
        end
        
        edges=linspace(min(pca_input_matrix(:,par_num)), max(pca_input_matrix(:,par_num)),num_histnorm_cols);
        counts=nan(length(logical_vals),length(edges));
        for i=1:length(logical_vals)
            n = histc(pca_input_matrix(sig_vector==logical_vals(i) & is_pitd==logical_vals(num_area) & monkey_vector,par_num),edges);
            counts(i,:)=n;%/sum(n);
        end
        
        h=bar(edges,counts','stacked');
        %colormap(cmap)
        this_cmap=colormap;
        this_cmap(1,:)=bar_colors(num_area, :);
        colormap(this_cmap); %freezeColors;
        for i=1:length(h)
            set(get(h(i),'Children'),'EdgeColor',bar_colors(num_area, :))
            %set(get(h(i),'Children'),'FaceColor',bar_colors(num_row, :))
        end
        %xlim([-.6 1])
        pbaspect([172 142 1])
        box off
        title(sprintf('%01.2g%% of %3.0f cells sign. in %s',100*sum(counts(1,:))/sum(counts(:)),sum(counts(:)),monkey_names{monkey_used}));
        xlabel(contrast_names{par_num})
        ylim([0 max(sum(counts,1))+1])
        line([0 0],ylim,'Color','black','LineStyle','--');

    end
end

save_figures_to_harbor=false;
if save_figures_to_harbor
    for i=1:length(f)
    filename=['Attn_stats_' area{i} '_' monkey_names{monkey_used}];
    %plot2svg(fullfile(figure_dir,[filename '.svg']),f(i));
    saveas(f(i),fullfile(figure_dir,filename),'fig');
    saveas(f(i),fullfile(figure_dir,filename),'epsc');
    end
end

%% Save figures

if save_figures_to_harbor
    addpath('../plot2svg/')
    figure_dir='/Freiwald/ppolosecki/harbor';
    fig_handles=get(0,'children');
    for i=1:length(fig_handles)
        filename=['population_stats_figure_' sprintf('%02.0f',fig_handles(i))];
        plot2svg(fullfile(figure_dir,[filename '.svg']),fig_handles(i));
        saveas(fig_handles(i),fullfile(figure_dir,filename),'fig');
    end
end

