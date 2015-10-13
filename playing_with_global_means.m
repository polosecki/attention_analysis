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
sem_tcs=squeeze(nanmean(dpca_input{align}.matrix,1))/sqrt(size(dpca_input{align}.matrix,1));
%%
xlims_used=[-1 0];
ylims_used=[0 5];
figure;
subplot(1,3,1)
hold all
plot(t,squeeze(mean(mean_tcs(:,1,[2 4]),3)))
plot(t,squeeze(mean(mean_tcs(:,3,[2 4]),3)))
xlim(xlims_used)
ylim(ylims_used)
pbaspect([1 1 1])

subplot(1,3,2)
hold all
plot(t,squeeze(mean(mean_tcs(:,[2 4],1),2)))
plot(t,squeeze(mean(mean_tcs(:,[2 4],3),2)))
xlim(xlims_used)
ylim(ylims_used)
pbaspect([1 1 1])

subplot(1,3,3)
hold all
plot(t,mean_tcs(:,1,1))
plot(t,mean_tcs(:,3,1))
plot(t,mean_tcs(:,1,3))
plot(t,mean_tcs(:,3,3))
xlim(xlims_used)
ylim(ylims_used)
pbaspect([1 1 1])
legend

