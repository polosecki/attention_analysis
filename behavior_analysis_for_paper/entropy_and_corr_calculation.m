close all
nbins=50;
n_perms=1e4;
save_fig=true;
nbins=120;
for j=1:2
    f(j)=figure;
end
% for j=1:2
%     g(j)=figure;
% end

for mm=1:2
    
    [RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',mm);
    %%
    
    %data_mat=[BRT_abs-RT_abs nBRT_abs-RT_abs];
    data_mat=[RT_abs-BRT_abs RT_abs-nBRT_abs];
    N_mat=zeros(nbins,nbins);
    edges=linspace(min(data_mat(:)),max(data_mat(:)),nbins)';
    
    v1=data_mat(:,1); v2=data_mat(:,2);
    true_corr=corr(data_mat(:,1),data_mat(:,2))
    edges=linspace(min(data_mat(:)),max(data_mat(:)),nbins)';
    
    for used_col=1:2
        if used_col==1
            idx=true(size(data_mat,1),1);
        elseif used_col==2;
            idx=false(size(data_mat,1),1);
        end
        n = histc([data_mat(idx,1); data_mat(~idx,2)],edges)/size(data_mat,1);
        idx2= n~=0;
        true_ent(used_col) = - sum(log2(n(idx2)).*n(idx2));
    end
    
    rand_ent=zeros(n_perms,2);
    rand_corr=zeros(n_perms,1);
    tic
    ncount=0;
    for i=1:n_perms
        %       idx=logical(randi([0,1],size(data_mat,1),1));
        %       rand_data=[v1(idx)  v2(idx); v2(~idx)  v1(~idx)];
        b=RT_abs(randperm(length(RT_abs)));
        %rand_data=[BRT_abs-b nBRT_abs-b];
        rand_data=[b-BRT_abs b-nBRT_abs];
        for j=1:2
            z{j} = histc(rand_data(:,j),edges)/size(data_mat,1);
            idx2= z{j}~=0;
            rand_ent(i,j) = - sum(log2(z{j}(idx2)).*z{j}(idx2));
        end
        [bb,bc]=hist3(rand_data,'edges',{edges edges});
        N_mat=(ncount*N_mat+bb*100/sum(bb(:)))/(ncount+1);
        nocunt=ncount+1;
        
        rand_corr(i) = corr(rand_data(:,1),rand_data(:,2));
        %rand_corr(i) = corr(rand_data(1:size(data_mat,1)),rand_data(size(data_mat,1)+1:end));
    end
    toc
    rand_ent=rand_ent(:);
    %close all
    %figure;
    %ah=axes;
    for j=1:2
        set(0,'currentFigure',f(j))
        ah(mm,j)=subplot(2,1,mm);
        if j==1
            [N,xc]=hist(ah(mm,j),rand_ent,30);
            bar(xc,N/sum(N)*100); box off
            set(get(ah(mm,j),'children'),'FaceColor',[.5 .5 .5])
            %xlim([true_ent(1)-.1 true_ent(2)+.1])
            xlim([3.8 5.1])
            colors_used=distinguishable_colors(2);
            for i=1:length(true_ent)
                line([true_ent(i) true_ent(i)],ylim,'color',colors_used(i,:))
            end
            xlabel('Entropy of RT distribution, relative to t_{PME\_OT} (bits)')
            ylabel('%')
            title(['M' num2str(mm)])
            lg=legend({'Shuffled labels','Cued Surface','Dist. Surface'});
            set(lg,'box','off')
        elseif j==2
            [N,xc]=hist(ah(mm,j),rand_corr,30);
            bar(xc,N/sum(N)*100); box off
            set(get(ah(mm,j),'children'),'FaceColor',[.5 .5 .5])
            colors_used=distinguishable_colors(2);
            for i=1:length(true_corr)
                line([true_corr(i) true_corr(i)],ylim,'color',colors_used(i,:))
            end
            %xlim([0.2 0.6])
            xlim([0 1])
            xlabel(['Correlation between RTs relative to ' char(10) 't_{PME\_OT} of cued and distractor surface']);            ylabel('%')
            title(['M' num2str(mm)])
            lg=legend({'Shuffled RT trial labels','Real RT trial labels'});
            set(lg,'box','off')
        end
        
    end
    %        set(0,'currentFigure',g(mm))
    options=[];
    %options.title= ['M' num2str(mm)];
    options.cmap=cmap_used;
    options.zlim=[0 .15];
    options.ylim=[0 11];
    options.smooth_N=true;
    options.xlim=options.ylim;
    options.xlabel= 'RT relative to distractor t_{PME\_OT} (s)';
    options.ylabel= 'RT relative to cued t_{PME\_OT} (s)';
    %options.diag_hist=diff(data_mat,1,2)/sqrt(2);
    g(mm)=plot_2d_backup(N_mat,bc,options);
end
if save_fig
    for j=1:length(f)
       % export_fig(['/Freiwald/ppolosecki/harbor/entropy_BRT_distribution' num2str(j) '.eps'],'-eps', '-transparent',f(j))
    end
        for mm=1:length(g)
        export_fig(['/Freiwald/ppolosecki/harbor/fake_correlation_BRT_distribution_monkey_' num2str(mm) '.eps'],'-eps', '-transparent',g(mm))
    end
end
