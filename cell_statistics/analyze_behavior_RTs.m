clc;clear all; close all;
load full_stats; addpath('..')
monkey={'Quincy','Michel'};
area={'PITd','LIP'};

%%First question: Is there a spatial bias a given monkey?

aa=2; %PITd or LIP sessions
mm=2; %Quincy or Michel
used_arrangement=2; %0deg or 45deg tilted stimulus array
RTs=cell(2,2,2);

for aa=1:2
    for mm=1:2
            use_indexes=false(size(full_stats(aa,mm).RT_data));
            sessions_included=[];
            for cell_no=1:length(full_stats(aa,mm).RT_data)
                if ~isempty(full_stats(aa,mm).RT_data(cell_no).session_index) & ~ismember(full_stats(aa,mm).RT_data(cell_no).session_index,sessions_included)
                    use_indexes(cell_no)=true;
                    sessions_included=[sessions_included;full_stats(aa,mm).RT_data(cell_no).session_index];
                end
            end
            
            for cell_no=1:length(full_stats(aa,mm).RT_data)
                if any(full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions==0) & use_indexes(cell_no)
                    original_coords=full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions;
                    [sorted_coords,sorted_idx]=sortrows(original_coords);
                    RTs{aa,mm,1}=[RTs{aa,mm,1} cellfun(@median,full_stats(aa,mm).RT_data(cell_no).RTs(sorted_idx))];
                    labels{1}=sorted_coords(1:4,2);
                elseif any(full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions==45) & use_indexes(cell_no)
                    original_coords=full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions;
                    [sorted_coords,sorted_idx]=sortrows(original_coords);
                    RTs{aa,mm,2}=[RTs{aa,mm,2} cellfun(@median,full_stats(aa,mm).RT_data(cell_no).RTs(sorted_idx))];
                    labels{2}=sorted_coords(1:4,2);
                end
            end
                    end
end
for used_arrangement=1:2
    for aa=1:2
        for mm=1:2
                        RTs{aa,mm,used_arrangement}=sort(RTs{aa,mm,used_arrangement},2);

        end
    end
end
%Two-way anova:
Nreps=size(RTs{aa,mm,used_arrangement},2)
bb=reshape(RTs{aa,mm,used_arrangement}',[4*Nreps 4]);
%rows of bb are brt directions (and replications), columns are surface
%positions
figure; imagesc(bb)
set(gca,'xtick',1:4);
set(gca,'ytick',(1:4)*Nreps-round(Nreps/2));
set(gca,'xticklabel',labels{used_arrangement});
set(gca,'yticklabel',labels{used_arrangement});
xlabel('Attended Surface')
ylabel('BRT direction')

figure;
load('CoolWarmMap')
meaned=reshape(mean(RTs{aa,mm,used_arrangement},2),[4 4]);
scale_semiwidth=2*std(meaned(:));
scale_center=mean(meaned(:))
imagesc(meaned)
colormap(cmap)
set(gca,'xtick',1:4);
set(gca,'ytick',1:4);
set(gca,'xticklabel',labels{used_arrangement});
set(gca,'yticklabel',labels{used_arrangement});
title(['RTs= ' num2str(mean(meaned(:))) ' +- ' num2str(std(meaned(:))/sqrt(numel(meaned)-1))])

xlabel('Attended Surface')
ylabel('BRT direction')

[p, table, stats]=anova2(bb,Nreps);
figure;
%c=multcompare(stats,'estimate','column');
c=multcompare(stats,'estimate','row');
%% OVERALL REACTION TIMES:
for aa=1:2
    for mm=1:2
    %RTs_collapsed{aa,mm}=vertcat(RTs{aa,mm,1}(:),RTs{aa,mm,2}(:));
    RTs_collapsed{aa,mm}=horzcat(mean(RTs{aa,mm,1},1),mean(RTs{aa,mm,2},1));

    end
end

cellfun(@mean,RTs_collapsed)
cellfun(@std,RTs_collapsed)
%% Behavioral question: are RTs affected by relative position of attended
% surface and BRT direction?
close all; 
RTs=cell(2,2,2);
aa=1;
%1-Find indexes of data from individual sessions (non-repeated):
use_indexes=false(size(full_stats(aa,mm).RT_data));
sessions_included=[];
for aa=1:2
    for mm=1:2
for cell_no=1:length(full_stats(aa,mm).RT_data)
    if ~isempty(full_stats(aa,mm).RT_data(cell_no).session_index) & ~ismember(full_stats(aa,mm).RT_data(cell_no).session_index,sessions_included)
        use_indexes(cell_no)=true;
        sessions_included=[sessions_included;full_stats(aa,mm).RT_data(cell_no).session_index];
    end
end

%2-Browse throught cells:
for cell_no=1:length(full_stats(aa,mm).RT_data)
    if any(full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions==0) & use_indexes(cell_no)
        original_coords=full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions;
        [sorted_coords,sorted_idx]=sortrows(original_coords);
        RTs{aa,mm,1}=[RTs{aa,mm,1} cellfun(@median,full_stats(aa,mm).RT_data(cell_no).RTs(sorted_idx))];
        labels{1}=sorted_coords(1:4,2);
    elseif any(full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions==45) & use_indexes(cell_no)
        original_coords=full_stats(aa,mm).RT_data(cell_no).conditions_abs_positions;
        [sorted_coords,sorted_idx]=sortrows(original_coords);
        RTs{aa,mm,2}=[RTs{aa,mm,2} cellfun(@median,full_stats(aa,mm).RT_data(cell_no).RTs(sorted_idx))];
        labels{2}=sorted_coords(1:4,2);
    end
end
        end
    end

%%
%Determine consistent, inconsistent and orthogonal conditions:
sorted_coords=allcomb(0:90:270,0:90:270);
is_consistent=(diff(sorted_coords,1,2)==0);
is_inconsistent=abs(diff(sorted_coords,1,2))==180;
is_orthogonal= ~is_consistent & ~is_inconsistent;
consistency_cases=[is_consistent is_inconsistent is_orthogonal];

RTs_by_consistency=cell(3,1);
consistency_labels=cell(size(RTs_by_consistency));
for cons_index=1:3
    for used_arrangement=1:2
        for aa=1:2
            garch{aa,used_arrangement}=RTs{aa,mm,used_arrangement}(consistency_cases(:,cons_index),:);
            %RTs_by_consistency{cons_index}=garch(:); clear garch;
        end
    end
    %keyboard
    RTs_by_consistency{cons_index}=vertcat(garch{1}(:),garch{2}(:),garch{3}(:),garch{4}(:)); clear garch
    consistency_labels{cons_index}=cons_index*ones(size(RTs_by_consistency{cons_index}));
end

%[p, table, stats]=anova1(vertcat(RTs_by_consistency{:}),vertcat(consistency_labels{:}));
[p, table, stats]=kruskalwallis(vertcat(RTs_by_consistency{1:2}),vertcat(consistency_labels{1:2}));

title(['Attention-Response consistency effect - Area: ' area{aa} ', Monkey: ' monkey{mm} ', Tilt: ' num2str(used_arrangement-1) ' ,p=' sprintf('%1.1g',p)])
%figure;
%c=multcompare(stats)

%%



