function [results]= make_GLM_and_contrasts_from_inst_firing(y,closest_surf_phi,surf_str,noise_type,mean_center,std_scale)

%% Define variables
do_anovan=0;
RF_surf_index=find(unique([surf_str.phi])'==closest_surf_phi);

locations=circshift(unique([surf_str.phi])',-RF_surf_index+1);%RF_centered; for computation purposes
%non_shifted_locations=unique([surf_str.phi])'; %For nomenclature purposes

locations=reshape(locations,2,2); %rows indicate physical conditions
%non_shifted_locations=reshape(non_shifted_locations,2,2);
phys_surf=2*ismember([surf_str.phi]',locations(1,:))-1;
phys_targ=2*ismember([surf_str.brt]',locations(1,:))-1;

attend1=zeros(size(phys_targ)); attend1(ismember([surf_str.phi]',locations(1,1)))=1;
attend1(ismember([surf_str.phi]',locations(1,2)))=-1;

attend2=zeros(size(phys_targ)); attend2(ismember([surf_str.phi]',locations(2,1)))=1;
attend2(ismember([surf_str.phi]',locations(2,2)))=-1;

%attend=[surf_str.phi];
%saccade=[surf_str.brt]';

saccad1=zeros(size(phys_targ)); saccad1(ismember([surf_str.brt]',locations(1,1)))=1;
saccad1(ismember([surf_str.brt]',locations(1,2)))=-1;

saccad2=zeros(size(phys_targ)); saccad2(ismember([surf_str.brt]',locations(2,1)))=1;
saccad2(ismember([surf_str.brt]',locations(2,2)))=-1;

%% Run ANOVAN
if do_anovan
    variable_names={'Surface Position','Target Position','Attention','Saccade'};
    conds={phys_surf,phys_targ,attend,saccade};
    
    % 'nested': A matrix M of 0's and 1's specifying the nesting relationships among the grouping variables.
    %M(i,j) is 1 if variable i is nested in variable j
    nested_matrix=zeros(length(variable_names));
    nested_matrix(3,1)=1; % var 3 (attention) is nested within var 1 (surf position)
    nested_matrix(4,2)=1; % var 4 (Saccade) is nested within var 2 (Target Position)
    
    %Statistical tests of all effects of interes are calculated below:
    close all;
    results.anovan.p=[];
    for i=1:size(y,2)
        disp(['ANOVAN ' num2str(i)])
        [p,table,stats,terms] = anovan(y(:,i),conds,'varnames',variable_names,'model','interaction','nested',nested_matrix);
        %stats.coeffs give you effect sizes, stats.coeffnames give you what is what
        close(1)
        results.anovan.p=[results.anovan.p p];
        results.anovan.table{1,i}=table;
        results.anovan.stats{1,i}=stats;
        results.anovan.stats{i}=stats;
        results.anovan.terms{i}=terms;
    end
end
%% RUN GLMs for effect sizes
%To obtain effect sizes am interested in:
%-Main effect of surface postion
%-Main effect of target postion
%-Main effect of saccade
%-Main effect of attention
%-Interaction attention / target
%-Interaction attention / saccade

%Interactions:
interaction_pairs={'phys_surf','phys_targ';
    'phys_surf', 'saccad1';
    'phys_surf', 'saccad2';
    'phys_targ', 'attend1';
    'phys_targ', 'attend2';
    'attend1', 'saccad1';
    'attend1', 'saccad2';
    'attend2', 'saccad1';
    'attend2', 'saccad2';};
interact_predictors=nan(length(phys_surf),size(interaction_pairs,1));
for c=1:size(interaction_pairs,1)
    eval_string=[interaction_pairs{c,1} '==max(' interaction_pairs{c,1} ') & ' interaction_pairs{c,2} '==max( ' interaction_pairs{c,2} ')'];
    temp=double(eval(eval_string));
    temp(temp==0)=-1;
    eval_string=[interaction_pairs{c,1} '==0 | ' interaction_pairs{c,2} '==0'];
    zero_indexes=eval(eval_string);
    temp(zero_indexes)=0;
    interact_predictors(:,c)=temp/range(temp);
end
X3=[phys_surf/2 phys_targ/2 attend1/2 attend2/2 saccad1/2 saccad2/2 ones(size(phys_surf)) interact_predictors];

cols_surf_only_in_RF=logical([0 0 1 0 0 1 1 0 0 0 0 0 0 1 0 0]);
cols_targ_only_in_RF=logical([0 0 0 1 1 0 1 0 0 0 0 0 0 0 1 0]);
all_cols=true(size(cols_surf_only_in_RF));
col_cells={cols_surf_only_in_RF,cols_targ_only_in_RF,all_cols};

trials_surf_only_in_RF= phys_surf==1  & phys_targ==-1;
trials_targ_only_in_RF= phys_surf==-1 & phys_targ==1;
all_trials=true(size(trials_surf_only_in_RF));

trials_used={trials_surf_only_in_RF,trials_targ_only_in_RF,all_trials};

% contrast.matrix={[1 0 0 0 0 0 0]; %Surface
%                  [0 1 0 0 0 0 0]; %Sacc target
%                  [0 0 1 0 0 0 0]; %Attention 1
%                  [0 0 0 1 0 0 0]; %Attention 2
%                  [0 0 0 0 1 0 0]; %Saccade 1
%                  [0 0 0 0 0 1 0]}; %Saccade 2

contrast.name={'Surface';
    'Target';
    'Attn(Surf in RF)';
    'Attn(Surf out RF)';
    'Saccade(Targ in RF)';
    'Saccade(Targ out RF)'
    'mean'};
%           keyboard
temp_names=cell(size(interaction_pairs,1),1);
for ii=1:size(interaction_pairs,1)
    temp_names{ii}=['int_' interaction_pairs{ii,1} '_' interaction_pairs{ii,2}];
end

contrast.name=[contrast.name; temp_names];

for used_mat=1:3
    X=X3(:,col_cells{used_mat});
    if used_mat==3
        [Q,R] = qr(X);
        D=diag(sign(diag(R)));
        orthogonal_X=(Q(:,1:size(X,2))*D).*repmat(sqrt(diag(X'*X)'),size(X,1),1);
        X(:,8:end)=orthogonal_X(:,8:end);
    end
    %crazy_mat3=inv(X'*X);
    results.GLM(used_mat).ces_std=zeros(size(X,2),size(y,2));
    results.GLM(used_mat).t=zeros(size(X,2),size(y,2));
    results.GLM(used_mat).Fsig=zeros(size(X,2),size(y,2));
    results.GLM(used_mat).ces=zeros(size(X,2),size(y,2));
    results.GLM(used_mat).dof=zeros(size(X,2),size(y,2));
    rx=rank(X);
    disp(['GLM(' num2str(used_mat) ')'])
    
    estpdisp='on';
    for i=1:size(y,2)
           disp(['GLM(' num2str(used_mat) ') ' num2str(i)])
        temp_y=y(:,i);
        tt=trials_used{used_mat} & ~isnan(temp_y);
        if rank(X(tt,:))<rx
            continue
        end   
        [b,dev,stats] = glmfit(X(tt,:),temp_y(tt),noise_type,'link','identity','constant','off','estdisp',estpdisp);
        results.GLM(used_mat).ces_std(:,i)=stats.se/std_scale;
        results.GLM(used_mat).t(:,i)=stats.t;
        results.GLM(used_mat).Fsig(:,i)=stats.p;
        results.GLM(used_mat).ces(:,i)=stats.beta/std_scale;
        results.GLM(used_mat).dof(:,i)=stats.dfe;
    end
    results.GLM(used_mat).X=X;
    results.GLM(used_mat).contrast.name=contrast.name(col_cells{used_mat});
    %results.GLM(used_mat).contrast.matrix=contrast.matrix;
end
%To check ortogonality of the matrix:
% sx=X3.^2;
% a=(sum(sx,1)).^(.5);
% normat=repmat(a,size(X3,1),1);
% X0=X3./normat;
% imagesc(X0'*X0)