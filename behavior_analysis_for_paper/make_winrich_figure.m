close all
%Centered on RT:

nbins=120;
cmap_used='bone';'hot';
save_fig=true;

for mm=1:2
[RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',mm);
data_mat=[RT_abs-BRT_abs RT_abs-nBRT_abs];
edges=linspace(min(data_mat(:)),max(data_mat(:)),nbins)';
[N,bc]=hist3(data_mat,'edges',{edges edges});
%[N,bc]=hist3([BRT_abs-RT_abs nBRT_abs-RT_abs],'edges',{edges edges});
N=N/sum(N(:))*100;
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
%f=plot_2d_distribution_with_marginals(N,bc,options);
f=plot_2d_backup(N,bc,options);

if save_fig
    export_fig(['/Freiwald/ppolosecki/harbor/triggered_joint_distribution_monkey_' num2str(mm) '.eps'],'-eps','-transparent',f)
end
end

%%
if false
%subplot(3,2,3)

    [RT_abs, BRT_abs, ~, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',1);
    edges=linspace(min(RT_abs-BRT_abs),max(RT_abs-BRT_abs),nbins)';
    [N,bc]=hist3([RT_abs-BRT_abs nBRT_abs-BRT_abs],'edges',{edges edges});
    N=N/sum(N(:))*100;
    
    options.ylabel='RT(s)';
    options.xlabel='Distractor BRT RT(s)';
    options.title='M1 - relative to cued BRT';
    
    f=plot_2d_distribution_with_marginals(N,bc,options)
    


%%
%subplot(3,2,5)
    
    [RT_abs, BRT_abs, ~, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',1);
    edges=linspace(min(RT_abs-nBRT_abs),max(RT_abs-nBRT_abs),nbins)';
    [N,bc]=hist3([RT_abs-nBRT_abs BRT_abs-nBRT_abs],'edges',{edges edges});
    N=N/sum(N(:))*100;
    
    options.ylabel='RT(s)';
    options.xlabel='Cued BRT RT(s)';
    options.title='M1- relative to distractor BRT';
    
    f=plot_2d_distribution_with_marginals(N,bc,options)
    
    %%
    %subplot(3,2,2)
    [RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',2);
    edges=linspace(min(BRT_abs-RT_abs),max(BRT_abs-RT_abs),nbins)';
    [N,bc]=hist3([BRT_abs-RT_abs nBRT_abs-RT_abs],'edges',{edges edges});
    N=N/sum(N(:))*100;
    
    options.xlabel='Dist Surface BRT(s)';
    options.ylabel='Cued Surface BRT(s)';
    options.title='M2 - relative to RT';
    
    f=plot_2d_distribution_with_marginals(N,bc,options)
    
    %%
    %subplot(3,2,4)
    [RT_abs, BRT_abs, ~, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',2);
    edges=linspace(min(RT_abs-BRT_abs),max(RT_abs-BRT_abs),nbins)';
    [N,bc]=hist3([RT_abs-BRT_abs nBRT_abs-BRT_abs],'edges',{edges edges});
    N=N/sum(N(:))*100;
    options.ylabel='RT(s)';
    options.xlabel='Distractor BRT RT(s)';
    options.title='M2 - relative to cued BRT';
    
    f=plot_2d_distribution_with_marginals(N,bc,options)
    
    %%
    %subplot(3,2,6)
    [RT_abs, BRT_abs, nRT_abs, nBRT_abs] = get_distribution_from_total_behav(total_behav,'none',2);
    edges=linspace(min(RT_abs-nBRT_abs),max(RT_abs-nBRT_abs),nbins)';
    [N,bc]=hist3([RT_abs-nBRT_abs BRT_abs-nBRT_abs],'edges',{edges edges});
    N=N/sum(N(:))*100;
    options.ylabel='RT(s)';
    options.xlabel='Cued BRT RT(s)';
    options.title='M2 - relative to distractor BRT';
    
    f=plot_2d_distribution_with_marginals(N,bc,options)
end
%%

