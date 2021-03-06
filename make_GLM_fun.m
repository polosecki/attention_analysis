function [results]=make_GLM_fun(cell_no,monkey,area,operation_mode,z_score_data,use_BRTs,noise_type)

%Inputs: cell_no: cell number in data base
%         monkey: monkey name
%         area: brain area
%         operation_mode: 'time_course' for full time analysis with plot
%                         'fixed_points' for just fixed set of time points (for histogram of population)
%                         'betas_for_pca' for high temporal resolution and
%                         no plot
%         z_score: set to 1 to normalize data
%         use_BRT: include BRT-centerered activity in results
%         noisetype: 'normal' or 'poisson'

%cell_no=28;%PITd; %28;% LIP;
if nargin<6
    use_BRTs=1;
end
if nargin<7
    noise_type='normal';
end

base_dir='/Freiwald/ppolosecki/lspace';
%if strcmp(monkey,'Michel') & strcmp(area,'LIP')
%    base_dir='/Freiwald/ppolosecki/lspace/sara';
%end
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
%monkey='Quincy';
%area='LIP';%'IP';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);

fixed_time_bins={[0.1 1.2 -0.5]; %stim-onset; mid-trial; pre-stim
    [-0.65 -.1]; %one seclong before saccade onset
    []}; %
fixed_time_semiwidths={[0.2 0.2 0.2]; %stim-onset; mid-trial; pre-stim
    [0.15 .1]; %one seclong before saccade onset
    []}; %

%% Load basic files

load(cell_file);
load(results_file)
load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.mat{end}));
ifname=fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.hdf{end});

fid = dhfun('open',ifname);
tmap = dh_get_trialmap_struct(fid);
dhfun('close',fid);

[surf_str] = trial_info(tmap,tds);
surf_str_extra = heiko_log_surf_params(tds{1},false); % used to be called p in michael code

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

if ~exist(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']),'file');
    error('Run the first atention analysis code to create the PSTHs')
else
    load(fullfile(base_dir,lower(monkey),cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']));
end


%% Run GLM and Tests:
%Commands of possible relevance:
%
%mvregress: multivariate general linear model
%manova1:multavariate analysis of variance
%fast_glmfit: fit glm model (FSFAST)
%fast_fratio: evaluate contrasts
%anovan: N-way anova

if use_BRTs
    num_matrices=length(grand_psth.matrix);
else
    num_matrices=2;
end

if z_score_data
    good_times={[-1 1],[-.8 0],[-.3 .3]};
    raw_grand=[];
    for i=1:num_matrices
        idx=(grand_psth.time_axis{i}>good_times{i}(1) & grand_psth.time_axis{i}<good_times{i}(2));
        raw_grand=[raw_grand grand_psth.matrix{i}(:,idx)];
    end
    conditions=allcomb(unique([surf_str.phi]'),unique([surf_str.brt]'));
    mean_mat=nan(size(conditions,1),size(raw_grand,2));
    zz=zeros(length(grand_psth.trials_used{1}),num_matrices);
    for col=1:num_matrices
        zz(:,col)=grand_psth.trials_used{col};
    end
    trials_used=prod(zz,2);
    for j=1:length(conditions)
        these_trials=trials_used & [surf_str.phi]'==conditions(j,1) & [surf_str.brt]'==conditions(j,2);
        mean_mat(j,:)=nanmean(raw_grand(these_trials,:));
    end
    
    mean_center=nanmean(mean_mat(:));
    std_scale=nanstd(mean_mat(:));
else
    mean_center=0;
    std_scale=1;
end

for mat_used=1:num_matrices
    
    %    tres = 1e9/30000; % Specified in nanoseconds, reciprocal of sample rate
    %
    %     % Define time:
    %     t = 0:tres*grand_psth.psthdec/1e9:(size(grand_psth.matrix{mat_used},2)-1)*tres*grand_psth.psthdec/1e9;
    %     if mat_used==1
    %         %Center on surface onset(stmulus onset+pause_duration):
    %         t_zero=-grand_psth.start_trig(mat_used).offset+unique(surf_str_extra.pausedur(grand_psth.trials_used{mat_used}));
    %     elseif mat_used==2
    %         %Center on saccade onset:
    %         t_zero=t(end)-grand_psth.end_trig(mat_used).offset;
    %     end;
    %     t=t-t_zero;
    
    t=grand_psth.time_axis{mat_used};
    %trial_count=sum(~isnan(grand_psth.matrix{mat_used}),1);
    %times_used=trial_count>15;
    RF_surf=nan(length(results_data),1);
    RF_surf(~cellfun(@isempty,{results_data.closest_surf_phi}))=[results_data.closest_surf_phi]';
    
    switch operation_mode
        case 'time_course'
            tbins_center=t(3:3:end); %50ms bins
            tbins_semi_width=0.15; % i.e., width=2*semi_width
        case 'fixed_points'
            tbins_center=fixed_time_bins{mat_used};
            tbins_semi_width=fixed_time_semiwidths{mat_used}; %variable width
        case 'betas_for_pca'
            tbins_center=t;
            tbins_semi_width=0;
        otherwise
            error('Set a valid Operation Mode')
            
    end
    y=nan(size(grand_psth.matrix{mat_used},1),length(tbins_center));
    switch operation_mode
        case 'fixed_points'
            for i=1:length(tbins_center)
                y(:,i)=nanmean(grand_psth.matrix{mat_used}(:,abs(t-tbins_center(i))<=tbins_semi_width(i)),2);
            end
        otherwise
            for i=1:length(tbins_center)
                y(:,i)=nanmean(grand_psth.matrix{mat_used}(:,abs(t-tbins_center(i))<=tbins_semi_width),2);
            end
    end    
    [temp]= make_GLM_and_contrasts_from_inst_firing(y,RF_surf(cell_no),surf_str,noise_type,mean_center,std_scale);
    results{mat_used}=temp; clear temp;
    results{mat_used}.time=tbins_center;
            y=(y-mean_center)/std_scale;   
    results{mat_used}.y=y;
    results{mat_used}.mean_activity=mean_center;
    results{mat_used}.std_activity=std_scale;
    results{mat_used}.noise_model=noise_type;
end
%% Make plots
switch operation_mode
    case 'time_course'
%         switch noise_type
%             case 'normal'
%                 contrasts_plotted={logical([1 1 0 0 0 0 1 0 0 1]);
%                     logical([0 0 0 0 1 0 1 1 0 0 1]);
%                     logical([0 0 0 0 0 0])};
%             case 'poisson'
                 contrasts_plotted={logical([0 0 0 0 0 0 0 0 0 0]);
                     logical([0 0 0 0 0 0 0 0 0]);
                     logical([1 1 1 0 1 0 0 1 1 0 1 0 1 0 0 0])};
%          end
        plot_GLM_contrasts(results,contrasts_plotted,use_BRTs,cell_str,cell_no);
    case 'betas_for_pca'
        %         contrasts_plotted={logical([0 0 0 0 0 0 0 0 0 0]);
        %                            logical([0 0 0 0 0 0 0 0 0]);
        %                            logical([1 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0])};
        
        contrasts_plotted={logical([1 1 0 0 0 0 1 0 0 1]);
            logical([0 0 0 0 1 0 1 1 0 0 1]);
            logical([0 0 0 0 0 0])};
        % plot_GLM_contrasts(results,contrasts_plotted,use_BRTs,cell_str,cell_no);
end


