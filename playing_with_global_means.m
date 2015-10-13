close all; clear all
monkey='both';%'both';%
area='LIP';

if strcmp(monkey,'Quincy') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=false;
    load dPCA/PITd_Quincy_dpca_input_with_mean.mat    
elseif strcmp(monkey,'Quincy') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=false;
    load dPCA/LIP_Quincy_dpca_input_with_mean.mat
elseif strcmp(monkey,'Michel') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=false;
    load dPCA/LIP_Michel_dpca_input_with_mean.mat
    base_dir='/Freiwald/ppolosecki/lspace/sara';
elseif strcmp(monkey,'Michel') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=false;
    load dPCA/PITd_Michel_dpca_input_with_mean.mat
end
if strcmp(monkey,'both')
    monkeys={'Quincy','Michel'};
    for mm=1:length(monkeys)
        load(['dPCA/' area '_' monkeys{mm} '_dpca_input_with_mean.mat'])
        if mm==1
            pt=dpca_input;
        end
    end
    for i=1:length(dpca_input)
        dpca_input{i}.matrix=[pt{i}.matrix;dpca_input{i}.matrix];
    end
end

align=2;
t=dpca_input{align}.time_axis;

mean_tcs=squeeze(nanmean(dpca_input{align}.matrix,1));
var_tcs=squeeze(nanvar(dpca_input{align}.matrix,1))/(size(dpca_input{align}.matrix,1));
%%

plot_params=containers.Map;
plot_params('t')=t;
plot_params('xlims')=[-.6 0.025];
plot_params('ylims')=[0 5];
plot_params('pbaspect')=[1 1 1];

colorlist=distinguishable_colors(4);

xlims_used=[-.8 0];
ylims_used=[0 5];
figure;


subplot(1,3,1)
plot_params('colors')=colorlist([1 3],:);

mean_mat=[squeeze(mean(mean_tcs(:,1,[2 4]),3))';
    squeeze(mean(mean_tcs(:,3,[2 4]),3))'];

strd_mat=[sqrt(squeeze(mean(var_tcs(:,1,[2 4]),3)))'/sqrt(2);
    sqrt(squeeze(mean(var_tcs(:,3,[2 4]),3)))'/sqrt(2)];
means_psth_plot(mean_mat,strd_mat,plot_params)



subplot(1,3,2)
plot_params('colors')=repmat(colorlist(4,:),2,1);
plot_params('linestyle')={'-',':'}

mean_mat=[squeeze(mean(mean_tcs(:,[2 4],1),2))';
    squeeze(mean(mean_tcs(:,[2 4],3),2))'];

strd_mat=[sqrt(squeeze(mean(var_tcs(:,[2 4],1),2)))'/sqrt(2);
    sqrt(squeeze(mean(var_tcs(:,[2 4],3),2)))'/sqrt(2)];
means_psth_plot(mean_mat,strd_mat,plot_params)

plot_params('colors')=repmat(colorlist([1 3],:),2,1);
plot_params('linestyle')={'-','-',':',':'};



subplot(1,3,3)
mean_mat=[mean_tcs(:,1,1)';
          mean_tcs(:,3,1)';
          mean_tcs(:,1,3)';
          mean_tcs(:,3,3)';
          ];
strd_mat=[sqrt(var_tcs(:,1,1))';
          sqrt(var_tcs(:,3,1))';
          sqrt(var_tcs(:,1,3))';
          sqrt(var_tcs(:,3,3))';
          ];
means_psth_plot(mean_mat,strd_mat,plot_params)

