close all; clear all
monkey='both';%'both';%
area='LIP';

%align=1; %1:align to surf onset; 2:align to saccade onset
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
base_dir='/Freiwald/ppolosecki/lspace';

cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

%% Load general files
addpath('..')

if ~strcmp(monkey,'both')
    load(cell_file);
    load(results_file)
    good_files=true(length(cell_str),1);
end

%Set good files and load z-cored, average PSTH of each condition, each cell
if strcmp(monkey,'Quincy') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=false;
    load ../dPCA/PITd_Quincy_dpca_input.mat    
elseif strcmp(monkey,'Quincy') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=false;
    load ../dPCA/LIP_Quincy_dpca_input.mat
elseif strcmp(monkey,'Michel') & strcmp(area,'LIP')
    good_files(bad_files(monkey,area))=false;
    load ../dPCA/LIP_Michel_dpca_input.mat
    base_dir='/Freiwald/ppolosecki/lspace/sara';
elseif strcmp(monkey,'Michel') & strcmp(area,'PITd')
    good_files(bad_files(monkey,area))=false;
    load ../dPCA/PITd_Michel_dpca_input.mat
end


if exist([area '_' monkey '_betas.mat'],'file')
    load([area '_' monkey '_betas.mat'])
end

if strcmp(monkey,'both')
    monkeys={'Quincy','Michel'};
    for mm=1:length(monkeys)
        load([area '_' monkeys{mm} '_betas.mat'])
        load(['../dPCA/' area '_' monkeys{mm} '_dpca_input.mat'])
        if mm==1
            tt=betas_population;
            pt=dpca_input;
        end
    end
    for i=1:length(betas_population)
        betas_population{i}.matrix=[tt{i}.matrix;betas_population{i}.matrix];
        dpca_input{i}.matrix=[pt{i}.matrix;dpca_input{i}.matrix];
    end
end





% To implement Bill Newsome approach (6.4-6.7 sections of supp material):
% %First: at the single cell level
% 1-Z-score each neurons response.
% 2-Make GLM fit with main exp variable beta values plus all pairwise interactions plus a time-dependent mean.
% %Second:
% 1-Run Population PCA to get rid of noise, reducing all responses to the first 12 PCs:
% PCAs are found for each time point and each condition.
% 2-For each tasks variable there is beta vector in neuron space.
% Each beta vector is denoised by projecting it to the first 12 PCs.
% For each deonised beta vector, a time is identified that maximizes its norm.
% Finally, the vectors are turned into an orthogonal basis using a QR decomposition.  (THIS STEP HAS DIFFERENT RESULTS DEPENDING ON THE  CHOICE OF ORDER FOR BETA VALUES)

% Useful functions/scipts:
% [linehandles]=plot_GLM_contrasts(results,contrasts_plotted)
% [results]=make_GLM_fun(cell_no,monkey,area,operation_mode)
% [results]= make_GLM_and_contrasts_from_inst_firing(y,closest_surf_phi,surf_str)
% scratch_attention_analysis_pp

%OBS: Z-SCORING BEFORE TAKING THE AVERAGE FOR EACH CONDITION CAN BE
%PROBLEMATIC, SINCE NON-BALANCED NUMBERS OF TRIALS MIGHT BIAS THING
%SLIGHTLY.

max_align=3;
use_single_effects_for_pca=true;
for align=1:max_align
    if use_single_effects_for_pca
        pca_input_matrix=reshape(dpca_input{align}.matrix(:,:,1:2,1:2),[size(dpca_input{align}.matrix,1) size(dpca_input{align}.matrix,2)*4])';
    else
        pca_input_matrix=reshape(dpca_input{align}.matrix,[size(dpca_input{align}.matrix,1) size(dpca_input{align}.matrix,2)*16])';
    end
    
% [coeff,score,latent,tsquared,explained,mu] = pca(X,Name,Value)
%coeff: says the pca coefficients
%score: each column of score represents one principal component
%latent: the variance explained by each component
%Xcentered = score*coeff'
%use: 'Centered',false to avoid the mean of each column (i.e., neuron) from being removed
%Num components: 'NumComponents',12 % to get 12 PCs only
[coeff,score,latent,~,explained,mu] = pca(pca_input_matrix,'Centered',false);

pca_elements(align).pca_input_matrix=pca_input_matrix;
pca_elements(align).coeff=coeff;
pca_elements(align).score=score;
pca_elements(align).explained=explained;
pca_elements(align).mu=mu;
end
%How to reconstruct the d_PCA matrix using nPCAs_used PCs:
nPCAs_used=12;
align=2;
X = pca_elements(align).score(:,1:nPCAs_used)*pca_elements(align).coeff(:,1:nPCAs_used)';

figure;
subplot(1,2,1)
imagesc(pca_elements(align).pca_input_matrix,[-2 5])
title('Original data')
subplot(1,2,2)
imagesc(X,[-2 5])
title(['De-noised data using ' num2str(nPCAs_used) ' PCs'])

%cell_no=29;
%%
t_boundaries={[-.5 1.5],[-1.5 0],[-.4 .4]}; % boundaries for each alignment
if ~exist('betas_population','var')
    for align=1:max_align
        betas_population{align}.matrix=zeros(sum(good_files),size(dpca_input{align}.matrix,2),4);
        betas_population{align}.time_axis=dpca_input{align}.time_axis;
    end
    
    for cell_no=1:length(cell_str)
        %cell_no=23%19
        if good_files(cell_no)
            i1=cumsum(good_files);
            disp(['Cell number: ' num2str(i1(cell_no)) ' out of ' num2str(max(i1))])
            [results]=make_GLM_fun(cell_no,monkey,area,'betas_for_pca',1);
            %close(gcf)
            for align=1:max_align
                t_start=t_boundaries{align}(1); %in secs, counting from t_zero (surf onset or saccade onset)
                t_end=t_boundaries{align}(2);
                % 3 is the analysis number for raw beta values.
                temp=results{align}.GLM(3).ces;
                t=results{align}.time;
                try_single_effects=true;
                if try_single_effects
                    %                temp(3,:)=results{align}.GLM(2).ces(8,:);
                    temp(3,:)=results{align}.GLM(1).ces(1,:);
                    %                temp(5,:)=2*results{align}.GLM(2).ces(10,:)-results{align}.GLM(2).ces(5,:);
                    temp(5,:)=results{align}.GLM(2).ces(2,:);
                end
                [~,index_start]=min(abs((t-t_start)));
                if cell_no==find(good_files,1)
                    %[~,index_end]=min(abs((t-t_end)));
                    index_duration(align)=round(diff(t_boundaries{align})/mode(diff(t)));
                    %else
                    %   index_end=index_start+index_duration(align);
                end
                index_end=index_start+index_duration(align);
                %[~,index_end]=min(abs((t-t_end)));
                temp=temp([3 5 1 2],index_start:index_end);%i. e.: attention,sacade,surface,target, in that order
                %figure; plot(t(index_start:index_end),temp');
                
                betas_population{align}.matrix(i1(cell_no),:,:)=temp';
            end
        end
    end
    save([area '_' monkey '_betas.mat'],'betas_population')
    
end
%% Denoise betas and calculate directions in PCA space:

nPCAs_used=6;
denoised_betas_population=betas_population;
for align=1:3
    cc=pca_elements(align).coeff(:,1:nPCAs_used); %PCA denoising matrix is cc*cc'
    figure;
    beta_directions{align}=zeros(size(betas_population{align}.matrix,1),size(betas_population{align}.matrix,3));
    for beta_index=1:size(betas_population{align}.matrix,3)
        beta_temp=squeeze(betas_population{align}.matrix(:,:,beta_index));
        denoised_temp=(cc*cc')*beta_temp;%beta_temp*(cc*cc');
        denoised_betas_population{align}.matrix(:,:,beta_index)=denoised_temp;
        beta_norm=sqrt(sum(denoised_temp.^2,1));
        [~,max_ind]=max(beta_norm);
%         if beta_index==2 && align==2%beta_index==3 && align==1%
%           %max_ind=round(length(beta_norm)/2);
%           max_ind=length(beta_norm);
%         end
        beta_directions{align}(:,beta_index)=denoised_temp(:,max_ind)/beta_norm(max_ind);
        cos_angle=beta_directions{align}(:,beta_index)'*(denoised_temp./repmat(beta_norm,size(denoised_temp,1),1));
        subplot(2,2,beta_index)
        plot(betas_population{align}.time_axis,100*(cos_angle.^2));
        ylim([0 100])
        %plot(100*sin(acos(cos_angle/max(cos_angle))))
        title(num2str(beta_index));
        xlabel('Time (s)')
    end
    % Orthogonalize the vector base using QR decomposition:
    
    [Q,R] = qr(beta_directions{align}(:,[2 1 3 4]));
    %[Q,R] = qr(beta_directions{align}(:,[1 2 3 4]));
    
    %[Q,R] = qr([ones(size(beta_directions{align},1),1) beta_directions{align}(:,[2 1 3 4])]);
    %[Q,R] = qr(beta_directions{align});
    D=diag(sign(diag(R))); % Fixes funny signs of the QR basis to be like in Gram-Schmidt
    beta_orthogonal_basis{align}=Q(:,1:size(betas_population{align}.matrix,3))*D;
    
    %beta_orthogonal_basis{align}(:,[1 2])=beta_orthogonal_basis{align}(:,[2 1]);
    %beta_orthogonal_basis{align}=Q(:,2:(size(betas_population{align},3)+1))*D(2:end,2:end);
end

%% Project Trajectories into the beta basis
align=2;
traces=zeros(size(dpca_input{align}.matrix,1),size(dpca_input{align}.matrix,2),2,2);
set(0,'DefaultAxesColorOrder',distinguishable_colors(20));
color_list=get(0,'DefaultAxesColorOrder');

plot_nums={{[1 2 3 4];[1 2 3 4]},{[1 3];[1 3]},{[1 2 3 4];[1 3]},{[1 3];[1 2 3 4]}};
%plot_nums={{[1 2 3 4];[1]},{[1 2 3 4];[2]},{[1 2 3 4];[3]},{[1 2 3 4];[4]}};

%f=figure; 
%set(f,'Position',get(0,'ScreenSize'));
save_figure=1;

for pp=1:length(plot_nums)
%subplot(2,2,pp); hold on
f=figure; hold on
%set(f,'Position',get(0,'ScreenSize'));

attn_conds=plot_nums{pp}{1};%[1 2 3 4];
sacc_conds=plot_nums{pp}{2};%[1 3];

axes_used=[1 2];%Axes: attn, sacc, surface, target
%axes_names={'Attention','Saccade','Surface presence','Sacc Target Presence'};
axes_names={'Sacade','Attention','Surface presence','Sacc Target Presence'};
%axes_names={'Attention','Surface presence','Sacade','Sacc Target Presence'};
cc=pca_elements(align).coeff(:,1:nPCAs_used);

legend_text={};
for aind=1:length(attn_conds)
for sind=1:length(sacc_conds)
    traces=dpca_input{align}.matrix;
    projected_denoised_trace=(beta_orthogonal_basis{align}')*(cc*cc')*traces(:,:,attn_conds(aind),sacc_conds(sind));
    projected_denoised_trace=projected_denoised_trace(axes_used,:);
    a=plot(projected_denoised_trace(1,:),projected_denoised_trace(2,:),'.-');
    %a=plot3(projected_denoised_trace(1,:),projected_denoised_trace(2,:),projected_denoised_trace(3,:),'.-');

    set(a,'color',color_list(length(sacc_conds)*(aind-1)+sind,:)) 
    legend_text{length(sacc_conds)*(aind-1)+sind}= ['attn=' num2str(attn_conds(aind)) '; sacc=' num2str(sacc_conds(sind))];
end
end
l=legend(legend_text,'location','NorthWest');%,'Orientation','horizontal')
set(l,'Box','off')
axis off
axis equal
xlabel(axes_names{axes_used(1)})
ylabel(axes_names{axes_used(2)})
title([monkey ' ' area ' target ' num2str(pp)] )



figure_dir='/Freiwald/ppolosecki/harbor';
i=0;
fig_filename=[area '_' monkey '_align_' num2str(align) '_' num2str(i)];

while exist(fullfile(figure_dir,[fig_filename '.svg']),'file')
    i=i+1;
    fig_filename=[fig_filename(1:end-1) num2str(i)];
end
if save_figure
    addpath('../plot2svg/')
    plot2svg(fullfile(figure_dir,[fig_filename '.svg']),f)
    saveas(f,fullfile(figure_dir,fig_filename),'fig');
end
end

%% Make figure for paper
align=2;
traces=zeros(size(dpca_input{align}.matrix,1),size(dpca_input{align}.matrix,2),2,2);
set(0,'DefaultAxesColorOrder',distinguishable_colors(20));
color_list=get(0,'DefaultAxesColorOrder');
line_styles={'-','--'};
marker_list={'.','none'};
marker_sizes=[27 4];
brightness=[1 0.7];%0.75]
plot_nums={{[1 3];[1 3]}};
%plot_nums={{[1 2 3 4];[1]},{[1 2 3 4];[2]},{[1 2 3 4];[3]},{[1 2 3 4];[4]}};

%f=figure; 
%set(f,'Position',get(0,'ScreenSize'));
save_figure=1;

for pp=1:length(plot_nums)
%subplot(2,2,pp); hold on
f=figure; hold on
set(f,'Position',get(0,'ScreenSize'));
%set (f, 'renderer', 'zbuffer')
attn_conds=plot_nums{pp}{1};%[1 2 3 4];
sacc_conds=plot_nums{pp}{2};%[1 3];

axes_used=[1 2];%Axes: attn, sacc, surface, target
%axes_names={'Attention','Saccade','Surface presence','Sacc Target Presence'};
axes_names={'Sacade','Attention','Surface presence','Sacc Target Presence'};
%axes_names={'Attention','Surface presence','Sacade','Sacc Target Presence'};
cc=pca_elements(align).coeff(:,1:nPCAs_used);

legend_text={};
for aind=1:length(attn_conds)
    for sind=1:length(sacc_conds)
        traces=dpca_input{align}.matrix;
        projected_denoised_trace=(beta_orthogonal_basis{align}')*(cc*cc')*traces(:,:,attn_conds(aind),sacc_conds(sind));
        projected_denoised_trace=projected_denoised_trace(axes_used,:);
        a=plot(projected_denoised_trace(1,:),projected_denoised_trace(2,:));
        %a=plot3(projected_denoised_trace(1,:),projected_denoised_trace(2,:),projected_denoised_trace(3,:),'.-');
        
        %set(a,'color',color_list(length(sacc_conds)*(aind-1)+sind,:))
        set(a,'color',color_list(1+length(sacc_conds)*(aind-1),:)*brightness(sind),'Marker',marker_list{sind},'MarkerFaceColor',color_list(1+length(sacc_conds)*(aind-1),:)*brightness(sind),'MarkerSize',marker_sizes(sind),'LineStyle',line_styles{sind},'LineWidth',2)
        legend_text{length(sacc_conds)*(aind-1)+sind}= ['attn=' num2str(attn_conds(aind)) '; sacc=' num2str(sacc_conds(sind))];
    end
end
l=legend(legend_text,'location','NorthWest');%,'Orientation','horizontal')
set(l,'Box','off')
axis off
axis equal
xlabel(axes_names{axes_used(1)})
ylabel(axes_names{axes_used(2)})
title([monkey ' ' area ' target ' num2str(pp)] )



figure_dir='/Freiwald/ppolosecki/harbor';
i=0;
fig_filename=[area '_' monkey '_align_' num2str(align) '_' num2str(i)];

while exist(fullfile(figure_dir,[fig_filename '.svg']),'file')
    i=i+1;
    fig_filename=[fig_filename(1:end-1) num2str(i)];
end
if save_figure
    addpath('../plot2svg/')
    plot2svg(fullfile(figure_dir,[fig_filename '.svg']),f)
    saveas(f,fullfile(figure_dir,fig_filename),'fig');
end
end