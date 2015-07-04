clear all; close all
monkey='Quincy';
area='PITd';
%cell_no=8;

overwrite_grand_psth=0;

plot_params(1).phi_used=[2 4];
plot_params(1).brt_used=1;
plot_params(1).align_used=[3 2];
plot_params(1).title='Saccade to RF, Surface away';

plot_params(2).phi_used=1;
plot_params(2).brt_used=1;
plot_params(2).align_used=[3 2];
plot_params(2).title='Saccade to RF, Attention in RF';

plot_params(3).phi_used=3;
plot_params(3).brt_used=1;
plot_params(3).align_used=[3 2];
plot_params(3).title='Saccade to RF, Attention away from RF';

plot_params(4).phi_used=1;
plot_params(4).brt_used=[2 4];
plot_params(4).align_used=[3 2];
plot_params(4).title='Attention in RF, Target away from RF';

% plot_params(5).phi_used=3;
% plot_params(5).brt_used=[2 4];
% plot_params(5).align_used=[3 2];
% plot_params(5).title='Attention away RF, Target away from RF';

save_figures=1;
send_text=true;

raster_fig_dir=fullfile('/Freiwald/ppolosecki/lspace/figures',[area '_' monkey]);
%% Load files
cell_file_dir='/Freiwald/ppolosecki/lspace/polo_preliminary/cell_file_manager';
cell_file=fullfile(cell_file_dir,[area '_' monkey '.mat']);
base_dir=fullfile('/Freiwald/ppolosecki/lspace/',lower(monkey));
results_file=fullfile(cell_file_dir,[area '_' monkey '_results.mat']);
load(cell_file)
load(results_file)

good_files=true(length(cell_str),1);
if strcmp(monkey,'Quincy') & strcmp(area,'PITd')
    good_files([4 8 22 23 24 32 33 42:46 52])=false;
elseif strcmp(monkey,'Quincy') & strcmp(area,'LIP')
    good_files([2 12:14 17 26])=false;
elseif strcmp(monkey,'Michel') & strcmp(area,'LIP')
    good_files([1 2 4 5 21 31 32 33])=0; % cell 1 removed beause it is a copy of cell 12
elseif strcmp(monkey,'Michel') & strcmp(area,'PITd')
    good_files([6 10:12 16 19 32 38 40 48])=0;
end
%%
for cell_no=1:length(cell_str)
    
    
    if good_files(cell_no)
        RF_pos=results_data(cell_no).RF_pos;
        
        if overwrite_grand_psth
            figure_dir=[];
            [closest_surf_phi] = attention_analysis_func_pp(cell_file,cell_no,base_dir,overwrite_grand_psth,figure_dir,RF_pos);
        end
        
        
        
        load(fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.mat{end}));
        ifname=fullfile(base_dir,cell_str(cell_no).dir,'proc',cell_str(cell_no).attention.hdf{end});
        
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
        
        if ~exist(fullfile(base_dir,cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']),'file');
            error('Run the first atention analysis code to create the PSTHs')
        else
            load(fullfile(base_dir,cell_str(cell_no).dir,'proc',['cell_' sprintf('%03.0f',cell_no) '_single_trial_attn_PSTH.mat']));
        end
        closest_surf_phi=results_data(cell_no).closest_surf_phi;
        RF_surf_index=find(unique([surf_str.phi])'==closest_surf_phi);
        locations=circshift(unique([surf_str.phi])',-RF_surf_index+1);%RF_centered; for computation purposes
        non_shifted_locations=unique([surf_str.phi])'; %For nomenclature purposes
        
        good_times={[-1 1],[-.8 0],[-.3 3]};
        raw_grand=[];
        for i=1:length(grand_psth)
            idx=(grand_psth.time_axis{i}>good_times{i}(1) & grand_psth.time_axis{i}<good_times{i}(2));
            raw_grand=[raw_grand grand_psth.matrix{i}(:,idx)];
        end
        conditions=allcomb(unique([surf_str.phi]'),unique([surf_str.brt]'));
        mean_mat=nan(size(conditions,1),size(raw_grand,2));
        zz=zeros(length(grand_psth.trials_used{1}),length(grand_psth));
        for col=1:length(grand_psth)
            zz(:,col)=grand_psth.trials_used{col};
        end
        trials_used=prod(zz,2);
        for j=1:length(conditions)
            these_trials=trials_used & [surf_str.phi]'==conditions(j,1) & [surf_str.brt]'==conditions(j,2);
            mean_mat(j,:)=nanmean(raw_grand(these_trials,:));
        end
        
        mean_center=nanmean(mean_mat(:));
        std_scale=nanstd(mean_mat(:));
        %%
        %align=2;
        close all
        
        
        
        RTs=[surf_str.rts]';
        [~,sorted_idxs]=sort(RTs);
        kern=ones(5,3); kern=kern/sum(kern(:));
        
        for plot_num=1:length(plot_params)
            f=figure;
            h=[];
            for align_num=1:length(plot_params(plot_num).align_used)
                this_plot_trials=grand_psth.trials_used{plot_params(plot_num).align_used(align_num)} & ismember([surf_str.phi]',locations(plot_params(plot_num).phi_used)) & ismember([surf_str.brt]',locations(plot_params(plot_num).brt_used));
                
                [~,idx]=min(abs(grand_psth.time_axis{plot_params(plot_num).align_used(align_num)}));
                plotted_mat=nanconv(grand_psth.matrix{plot_params(plot_num).align_used(align_num)}(sorted_idxs(this_plot_trials(sorted_idxs)),:),kern);
                
                if align_num==length(plot_params(plot_num).align_used)
                    
                end
                subplot(2,length(plot_params(plot_num).align_used),align_num)
                
                %imagesc(grand_psth.matrix{align}(this_plot_trials,:))
                %imagesc(plotted_mat)
                h(align_num)=image(nantowhite(plotted_mat,[mean_center-4*std_scale mean_center+4*std_scale],colormap('Winter')));
                set(gca,'XTick',[])
                line([idx idx],ylim,'Color', 'k')
                %axis off
                
                subplot(2,length(plot_params(plot_num).align_used),length(plot_params(plot_num).align_used)+align_num)
                plot(grand_psth.time_axis{plot_params(plot_num).align_used(align_num)},nanmean(plotted_mat))
                line([0 0],ylim,'Color', 'k')
                axis tight
                
                
            end
            
            [axh,labelh]=suplabel(plot_params(plot_num).title,'t');
            set(labelh,'FontSize',15);
            [~,temp,~]=fileparts(cell_str(cell_no).attention.hdf{end});
            [axh,labelh]=suplabel([temp '   cell_no: ' sprintf('%02.0f',cell_no)] ,'x');
            set(labelh,'Interpreter','none');
            %         for i=length(h)
            %             parent=get(h(i),'parent')
            %             set(parent,'Visible','off')
            %         end
            
            if save_figures
                filename=['cell_' sprintf('%03.0f',cell_no) '_raster_' sprintf('%01.0f',plot_num)];
                plot2svg(fullfile(raster_fig_dir,[filename '.svg']),f)
                saveas(f,fullfile(raster_fig_dir,filename),'png');
                saveas(f,fullfile(raster_fig_dir,filename),'fig');
            end
        end
    end
end
if send_text
    sendmail('9178422510@txt.att.net', 'Raster plots report', 'Loop completed');
end



%%

%%% Separate trials by BRt happening before or after the distractor BRT. 
% Perhaps restort to find_BRT_times_hieko_PP for this
%cued_length=cellfun(@length,bit_durs.cued);
%attended_brt_before_unattended_brt=(cellfun(@length,surf_str_extra.ncuedbitdurs)-cellfun(@length,surf_str_extra.cuedbitdurs))>0;


