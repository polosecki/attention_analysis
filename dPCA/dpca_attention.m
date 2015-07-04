clc; close all; clear all


monkey='Quincy';
area='PITd';

max_align=3; %set to 3 if BRT aligned stuff exists %1:align to surf onset; 2:align to saccade onset
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
base_dir='/Freiwald/ppolosecki/lspace';

cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

%% Load general files
addpath('..')
load(cell_file);
load(results_file)
good_files=true(length(cell_str),1);
if strcmp(monkey,'Quincy') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=false;
elseif strcmp(monkey,'Quincy') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=false;
elseif strcmp(monkey,'Michel') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=0; % cell 1 removed beause it is a copy of cell 12
elseif strcmp(monkey,'Michel') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=0;
end

%% Load cell-specific files for each cell
for cell_no=1:length(cell_str)
    if good_files(cell_no)
        load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.mat{end}));
        ifname=fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.hdf{end});
        
        fid = dhfun('open',ifname);
        tmap = dh_get_trialmap_struct(fid);
        dhfun('close',fid);
        
        [surf_str] = trial_info(tmap,tds);
        surf_str_extra = heiko_log_surf_params(tds{1},false); % used to be called p iichael code
        
        %---
        %Sanity checks:
        if any(strcmp({surf_str.out},'Success')'~=(tmap.oc==7))
            error('Tmap and surface info are not consistent')
        end
        test1=surf_str_extra.cuedsurfnames;
        test2=[surf_str.name];
        if ~all(strcmp(test1,test2))
            error('The surf_str surface info structures are not consistent')
        end
        %-----
        %Load grand psth matrix
        if ~exist(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']),'file');
            error('Run the first atention analysis code to create the PSTHs')
        else
            load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']));
        end
        
        if ~isequal(grand_psth.trials_used{1},grand_psth.trials_used{2})
            error('Your trials_used vectors are not equal. Re-run the creation of grand_psth')
        end
        
        
        %Select time window:
        t_boundaries={[-.5 1.5],[-1.5 0],[-.4 .4]}; % boundaries for each alignment
        
        for align=1:max_align
            t_start=t_boundaries{align}(1); %in secs, counting from t_zero (surf onset or saccade onset)
            t_end=t_boundaries{align}(2);  %in secs, counting from t_zero (surf onset or saccade onset)
            
            t=grand_psth.time_axis{align};
            [~,index_start]=min(abs((t-t_start)));
            if cell_no==find(good_files,1)                
            %[~,index_end]=min(abs((t-t_end)));
            %index_duration(align)=index_end-index_start;
            index_duration(align)=round(diff(t_boundaries{align})/mode(diff(t)));

            %else
             %   index_end=index_start+index_duration(align);
            end
               index_end=index_start+index_duration(align);
            this_trials=~isnan(grand_psth.matrix{align}(:,index_start:index_end));
            
            
            % Take the averages for each exp conditions
            phi_vect=[surf_str.phi]';
            brt_vect=[surf_str.brt]';
            closest_surf_RF=results_data(cell_no).closest_surf_phi; %get coords of surf closes to RF
            if any(phi_vect==45) %subtract 45 deg offsets on some cells
                phi_vect=phi_vect-45;
                brt_vect=brt_vect-45;
                closest_surf_RF=closest_surf_RF-45;
            end
            
            %Shift experimental conditions so that they are called in
            %RF-centered order
            min_index=find(unique(phi_vect)==closest_surf_RF);
            conditions=allcomb(unique(phi_vect),unique(brt_vect));
            conditions_shifted=allcomb(circshift(unique(phi_vect),-min_index+1),circshift(unique(brt_vect),-min_index+1));
            
            t_used=t(index_start:index_end);
            mean_mat=nan(size(conditions,1),length(t_used));
            for j=1:length(conditions)
                %Notice how the use of conditions_shifted selects the
                %trials in RF-specific coords
                this_trials=grand_psth.trials_used{align} & phi_vect==conditions_shifted(j,1) & brt_vect==conditions_shifted(j,2);
                mean_mat(j,:)=nanmean(grand_psth.matrix{align}(this_trials,index_start:index_end));
            end
            
            %Z-score the mean matrix
            %flattened_grand_mean=grand_psth.matrix{align}(this_trials,index_start:index_end);
            mm=nanmean(mean_mat(:));
            mstd=nanstd(mean_mat(:));
            mean_mat=(mean_mat-mm)/mstd;
            % Make multidimensional matrix for dPCA
            i1=cumsum(good_files); %index of cells in the good_file space
            phi_vals=unique(phi_vect);
            brt_vals=unique(brt_vect);
            if align==1 %make the 4-D matrix if it doesnt exist; note length of t_used varies fro each align value
                if ~exist('temp1','var')
                    temp1=nan(sum(good_files),length(t_used),length(unique(phi_vect)),length(unique(brt_vect)));
                end
            elseif align==2
                if ~exist('temp2','var')
                    temp2=nan(sum(good_files),length(t_used),length(unique(phi_vect)),length(unique(brt_vect)));
                end
            elseif align==3
                if ~exist('temp3','var')
                    temp3=nan(sum(good_files),length(t_used),length(unique(phi_vect)),length(unique(brt_vect)));
                end
            end
            
            dpca_input{align}.time_axis=t_used;
            for phi_ind=1:length(phi_vals)
                for brt_ind=1:length(brt_vals)
                    row_used=find(conditions(:,1)==phi_vals(phi_ind) & conditions(:,2)==brt_vals(brt_ind));
                    if align==1
                        temp1(i1(cell_no),:,phi_ind,brt_ind)=mean_mat(row_used,:);                        
                    elseif align==2
                        temp2(i1(cell_no),:,phi_ind,brt_ind)=mean_mat(row_used,:);
                    elseif align==3
                        temp3(i1(cell_no),:,phi_ind,brt_ind)=mean_mat(row_used,:);
                    end
                end
            end
            
            
        end
    end
end

dpca_input{1}.matrix=temp1;
dpca_input{2}.matrix=temp2;
dpca_input{3}.matrix=temp3;
save([area '_' monkey '_dpca_input'],'dpca_input');

%% Usa dPCA to find dimension that maximally explains variance accross exp conditions

clear comp; clear temp
n_components_used=5;

align=2;
temp=dpca_input{align}.matrix;

%Make dPCA base
W=dpca(temp,n_components_used,[],[]);

%Make projecion of each experimental conditions on each dPC
for stim2=1:4 
for stim1=1:4 
kk=W'*temp(:,:,stim1,stim2);
for comp_no=1:n_components_used
comp{comp_no}(stim1,stim2,:)=kk(comp_no,:);
end
end
end

%Make plots, one per dpc, with the projection of each condition plotted, to
%reveal wich exp variable each one is most associated with:
titled={'attn','sacc'}
for comp_no=1:n_components_used
    figure
    for align=1:2
    subplot(1,2,align)
    plot(squeeze(mean(comp{comp_no},3-align))')
    legend('1','2','3','4')
    title(['component: ' num2str(comp_no) ' ' titled{align}])
    end
end


%It seems dPC1 and dPC2 reflect attention and saccade direction
%respectively:
colors_used=distinguishable_colors(16);%varycolor(16);
comps_plotted=[1 2 3];
figure;
legend_vals={};
for stim1=1:4;%2:3
    for stim2=1:4
        if length(comps_plotted)==3
            h=plot3(squeeze(comp{comps_plotted(1)}(stim1,stim2,:)),squeeze(comp{comps_plotted(2)}(stim1,stim2,:)),squeeze(comp{comps_plotted(3)}(stim1,stim2,:)),'color',colors_used(stim1+4*(stim2-1),:));
            zlabel(['dPC: ' num2str(comps_plotted(3))])

        elseif length(comps_plotted)==2
            h=plot(squeeze(comp{comps_plotted(1)}(stim1,stim2,:)),squeeze(comp{comps_plotted(2)}(stim1,stim2,:)),'color',colors_used(stim1+4*(stim2-1),:));
        end
        set(h,'Marker','.')
        legend_vals{end+1}=['attn: ' num2str(stim1) ' sacc: ' num2str(stim2)];clc
        hold on
    end
end
grid on
xlabel(['dPC: ' num2str(comps_plotted(1))])
ylabel(['dPC: ' num2str(comps_plotted(2))])


legend(legend_vals,'location','EastOutside')
title([monkey ' ' area])






