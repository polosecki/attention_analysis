function [results]= make_GLM_and_contrasts_from_inst_firing(y,closest_surf_phi,surf_str)

%% Define variables
do_anovan=0;
RF_surf_index=find(unique([surf_str.phi])'==closest_surf_phi);

locations=circshift(unique([surf_str.phi])',-RF_surf_index+1);%RF_centered; for computation purposes
non_shifted_locations=unique([surf_str.phi])'; %For nomenclature purposes

locations=reshape(locations,2,2); %rows indicate physical conditions
non_shifted_locations=reshape(non_shifted_locations,2,2);
phys_surf=2*ismember([surf_str.phi]',locations(1,:))-1;
phys_targ=2*ismember([surf_str.brt]',locations(1,:))-1;

attend1=zeros(size(phys_targ)); attend1(ismember([surf_str.phi]',locations(1,1)))=1;
attend1(ismember([surf_str.phi]',locations(1,2)))=-1;

attend2=zeros(size(phys_targ)); attend2(ismember([surf_str.phi]',locations(2,1)))=1;
attend2(ismember([surf_str.phi]',locations(2,2)))=-1;

attend=[surf_str.phi];
saccade=[surf_str.brt]';

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


%Special vaiables for main effects and firstv interaction (attention / target position)
for i=1:2
    for j=1:2
att1_targ{i,j}=ismember([surf_str.phi]',locations(1,i)) &  ismember([surf_str.brt]',locations(j,:));
att2_targ{i,j}=ismember([surf_str.phi]',locations(2,i)) &  ismember([surf_str.brt]',locations(j,:));
    end
end
%att1_targ{end}=[];
%att2_targ{end}=[];ra
%For main

temp1=att1_targ';
temp2=att2_targ';
X1=double([[temp1{:}] [temp2{:}] saccad1 saccad2]);

%i.e., [att1_11 att1_12 att1_21 att1_22 att2_11 att2_12 att2_21 att2_22 sacc1 sacc2];
%phys surf = sum([temp1{:}],2) => sum([temp1{:}],2) - sum([temp2{:}],2) gives the effect of surface postion
%targ effect: [(att1_11 -att1_12 att1_21 -att1_22) (att2_11 -att2_12 att2_21 -att2_22)]
%atten1 effect: [(att1_11 att1_12 -att1_21 -att1_22)];
%atten2 effect: (att2_11 att2_12 -att2_21 -att2_22);
%atten1_taget_inter: [(att1_11 -att1_12 -att1_21 att1_22)]; 
%atten1_taget_inter: similar as above

contrast.matrix={[1 1 1 1 -1 -1 -1 -1 0 0]/4; %main effect of surface (possibly than this is better factoring atention out)
                 [1 -1 1 -1 1 -1 1 -1 0 0]/4; %main effect of target (would be good to have attention factore)
                 %[1 1 -1 -1 0  0 0  0 0 0]/2; %main effect of attnetion1
                 [1 0 -1 0 0  0 0  0 0 0]; % attention1 (w saccade targets in RF)
                 [0 1  0 -1 0  0 0  0 0 0]; % attention1 (w/no the saccade targets in RF)
                 [0 0  0 0  1 0 -1 0 0 0]; %attention2 (w saccade target in RF)
                 [0 0  0 0  0 1 0 -1 0 0]; %attention2 (w/no saccade target in RF) 
                 %[1 1 -1 -1 0  0 0  0 0 0;  %main effect of attnetion global
                 % 0 0  0 0  1 1 -1 -1 0 0]/2;
                 [1 -1 -1 1 0 0  0  0 0 0]/2;%att1/target interaction
                 [0  0 0  0 1 -1 1 -1 0 0]/2;%att2/target interaction
                 [0  0 0  0 0 0 1 1 0 0]/2; %baseline!! (i.e., nothing in RF)
                 [1 -1 1 -1 -1 1 -1 1 0 0]/4 %surface/target interaction
                 };
             
contrast.name={'Surface(main)';
               'Target (main)';
               %'Main effect of Attn(Stim in RF)';
               'Attn(Surf in RF)(Targ in RF)';
               'Attn(Surf in RF)(Targ out RF)';
               'Attn(Surf out RF)(Targ in RF)';
               'Attn(Surf out RF)(Targ out RF)';
               %'Main effect of attention, global';
               'Attn(Surf in RF)/Targ inter';
               'Attn(Surf out RF)/Targ inter'
               'Baseline'
               'Surf/Targ inter'};
 crazy_mat1=inv(X1'*X1);
 results.GLM(1).ces_std=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces_covar=cell(length(contrast.matrix),size(y,2));
 results.GLM(1).F=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).Fsig=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces_vector=cell(length(contrast.matrix),size(y,2));
 results.GLM(1).dof=zeros(length(contrast.matrix),size(y,2));
 for i=1:size(y,2)
     disp(['GLM(1) ' num2str(i)])
     temp_y=y(:,i);
     if rank(X1(~isnan(temp_y),:)'*X1(~isnan(temp_y),:))<rank(X1)
         continue
     end
     [beta, rvar, ~, ~] = fast_glmfit(temp_y(~isnan(temp_y)),X1(~isnan(temp_y),:));
     
     for c=1:length(contrast.matrix)
         [F, Fsig, ces, edof] = fast_fratio(beta,X1(~isnan(temp_y),:),rvar,contrast.matrix{c});
         ces_covar=rvar*contrast.matrix{c}*crazy_mat1*contrast.matrix{c}';
         ces_std=sqrt(mean(diag(ces_covar)))/sqrt(length(diag(ces_covar)));
         results.GLM(1).ces_std(c,i)=ces_std;
         results.GLM(1).ces_covar{c,i}=ces_covar;
         results.GLM(1).F(c,i)=F;
         results.GLM(1).Fsig(c,i)=Fsig;
         results.GLM(1).ces(c,i)=mean(ces);
         results.GLM(1).ces_vector{c,i}=ces;
         results.GLM(1).dof(c,i)=edof;
     end
     
     
     results.GLM(1).beta(:,i)=beta;
     results.GLM(1).rvar(1,i)=rvar;
     
 end
 results.GLM(1).X=X1;
 results.GLM(1).contrast.name=contrast.name;
 results.GLM(1).contrast.matrix=contrast.matrix;
 % Special variables for main effect of attn1 / saccade1 interaction

%i.e., [target att1_11 att1_12 att1_21 att2_22 att2 sacc2];

for i=1:2
    for j=1:2
att1_sacc1{i,j}=ismember([surf_str.phi]',locations(1,i)) &  ismember([surf_str.brt]',locations(1,j));
att2_sacc1{i,j}=ismember([surf_str.phi]',locations(2,i)) &  ismember([surf_str.brt]',locations(1,j));
    end
end
att1_bis=zeros(size(phys_surf)); %attention in RF, with no saccade targets in RF
att1_bis(ismember([surf_str.phi]',locations(1,1)) & ismember([surf_str.brt]',locations(2,:)))=1;
att1_bis(ismember([surf_str.phi]',locations(1,2)) & ismember([surf_str.brt]',locations(2,:)))=-1;

att2_bis=zeros(size(phys_surf)); %attention out of RF, with no saccade targets in RF
att2_bis(ismember([surf_str.phi]',locations(2,1)) & ismember([surf_str.brt]',locations(2,:)))=1;
att2_bis(ismember([surf_str.phi]',locations(2,2)) & ismember([surf_str.brt]',locations(2,:)))=-1;

temp=att1_sacc1';
temp2=att2_sacc1';


%TODO: an interaction that captures the interesting effect in LIP cell 11
%in Quincy
%X2_shitty=[phys_surf phys_targ [temp{:}] att1_bis attend2 saccad2 ones(size(phys_targ))]; %
X2=[phys_targ [temp{:}] [temp2{:}] att1_bis att2_bis saccad2];


%main effect of surface: no
 contrast2.matrix={2*[1 0 0 0   0 0 0 0 0 0 0 0]; %main effect of target
                   [0 1 1 -1 -1 0 0 0 0 0 0 0]/4 + ... %main effect of attention1:
                   2*[0 0 0 0 0 0 0  0  0 1 0 0]/2;
                     [0 0 0 0 0 1 1 -1 -1 0 0 0]/4 + ... %main effect of attention2:
                   2*[0 0 0 0 0 0 0  0  0 0 1 0]/2; 
                   2*[0 0 0 0 0 0 0  0  0 0 1 0]; % single effect of attention2 (with NOTHING in RF)
                     [0 0 0 0 0 1 -1 1 -1 0 0 0]/2;  % single effect of saccade 1* (no surf in RF)
                   2*[0 0  0 0  0 0  0 0  0 0 0 1];  % main effect of saccade 2
                     [0 1 -1 -1 1 0  0 0  0 0 0 0]/2; %attn1/sacc1 interaction*
                   2*[0 0  0  0 0 0  0 0  0 1 0 0]; %Atten1 (no target in RF)*
                     [0 1  1 -1 -1 0 0 0  0 0 0 0]/2; %Atten1 (target in RF)
                     [0 1 -1 1 -1 0 0  0 0 0 0 0]/4 + ...;%main effect of saccade1
                     [0 0  0 0  0 1 -1 1 -1 0 0 0]/4;
                     [0 1 -1 1 -1 -1 1 -1 1 0 0 0]/4; %Saccad1*/ Surface interaction
                     };

% contrast2.matrix={2*[1 0 0 0  0  0 0 0 0 0]; %main effect of surface 
%                   2*[0 1 0 0  0  0 0 0 0 0]; %main effect of target 
%                   [0 0 1 1 -1 -1 0 0 0 0]/4 + ... %main effect of attnetion1
%                   2*[0 0 0 0  0  0 1 0 0 0]/2;
%                   2*[0 0 0 0  0  0 0 1 0 0];%main effect of attention2
%                   [0 0 1 -1 1 -1 0 0 0 0]/2;%main effect of saccade1 *
%                   2*[0 0 0  0 0  0 0 0 1 0];%main effect of saccade2
%                   [0 0 1 -1 -1 1 0 0 0 0]/2; %attn1/sacc1 interaction*
%                   2*[0 0 0 0  0  0 1 0 0 0]; %Atten1 (no target in RF)*
%                   [0 0 1 1 -1 -1 0 0 0 0]/2; %Atten1 (target in RF)
%                   };


contrast2.name={'Target (main)';
               'Attn(Surf in RF) (main)';
               'Attn(Surf out RF) (main)';
               'Attn(Surf out RF) (targ out RF)';
               'Saccade(targ in RF) (Surf out RF)';
               'Saccade(targ out RF) (main)';
               'Attn(Surf in RF)/Sacc(Targ in RF) inter';
               'Attn(Surf in RF)(Targ out RF)';
               'Attn(Surf in RF)(Targ in RF)'
               'Saccade(targ in RF) (main)';
               'Saccade(targ in RF)/ Surf inter'};
           
% contrast2.name={'Surface (main)';
%                'Target (main)';
%                'Attn(Surf in RF) (main)';
%                'Attn(Surf out RF) (main)';
%                'Saccade(targ in RF) (main)';
%                'Saccade(targ out RF) (main)';
%                'Attn(Surf in RF)/Sacc(Targ in RF) inter';
%                'Attn(Surf in RF)(Targ out RF)'
%                'Attn(Surf in RF)(Targ in RF)'};
           

 results.GLM(2).ces_std=zeros(length(contrast2.matrix),size(y,2));
 results.GLM(2).ces_covar=cell(length(contrast2.matrix),size(y,2));
 results.GLM(2).F=zeros(length(contrast2.matrix),size(y,2));
 results.GLM(2).Fsig=zeros(length(contrast2.matrix),size(y,2));
 results.GLM(2).ces=zeros(length(contrast2.matrix),size(y,2));
 results.GLM(2).ces_vector=cell(length(contrast2.matrix),size(y,2));
 results.GLM(2).dof=zeros(length(contrast2.matrix),size(y,2));
 
% X2=[phys_surf phys_targ [temp{:}] att1_bis attend2 saccad2 ones(size(phys_targ))]; %
 crazy_mat2=inv(X2'*X2);
 for i=1:size(y,2)
     disp(['GLM(2) ' num2str(i)])
     temp_y=y(:,i);
     if rank(X2(~isnan(temp_y),:)'*X2(~isnan(temp_y),:))<rank(X2)
         continue
     end
     [beta2, rvar, ~, ~] = fast_glmfit(temp_y(~isnan(temp_y)),X2(~isnan(temp_y),:));
     
     for c=1:length(contrast2.matrix)
         [F, Fsig, ces, edof] = fast_fratio(beta2,X2(~isnan(temp_y),:),rvar,contrast2.matrix{c});
         ces_covar=rvar*contrast2.matrix{c}*crazy_mat2*contrast2.matrix{c}';
         ces_std=sqrt(mean(diag(ces_covar)))/sqrt(length(diag(ces_covar)));
         results.GLM(2).ces_std(c,i)=ces_std;
         results.GLM(2).ces_covar{c,i}=ces_covar;
         results.GLM(2).F(c,i)=F;
         results.GLM(2).Fsig(c,i)=Fsig;
         results.GLM(2).ces(c,i)=mean(ces);
         results.GLM(2).ces_vector{c,i}=ces;
         results.GLM(2).dof(c,i)=edof;
     end
     
     results.GLM(2).beta(:,i)=beta2;
     results.GLM(2).rvar(1,i)=rvar;
     
 end
 results.GLM(2).X=X2;
 results.GLM(2).contrast.name=contrast2.name;
 results.GLM(2).contrast.matrix=contrast2.matrix;
 
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
%keyboard
[Q,R] = qr(X3);
D=diag(sign(diag(R)));
orthogonal_X=(Q(:,1:size(X3,2))*D).*repmat(sqrt(diag(X3'*X3)'),size(X3,1),1);
X3(:,8:end)=orthogonal_X(:,8:end);

% contrast.matrix={[1 0 0 0 0 0 0]; %Surface
%                  [0 1 0 0 0 0 0]; %Sacc target
%                  [0 0 1 0 0 0 0]; %Attention 1
%                  [0 0 0 1 0 0 0]; %Attention 2
%                  [0 0 0 0 1 0 0]; %Saccade 1
%                  [0 0 0 0 0 1 0]}; %Saccade 2
contrast.matrix={[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Surface
                 [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Sacc target
                 [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Attention 1
                 [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]; %Attention 2
                 [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]; %Saccade 1
                 [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]; %Saccade 2
                 [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]; %Atten 1 - Saccade 1 inter
                 [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
                 [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
                 [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
                 [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];
                 [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];
                 [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]};
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
 crazy_mat3=inv(X3'*X3);
 results.GLM(3).ces_std=zeros(length(contrast.matrix),size(y,2));
 results.GLM(3).ces_covar=cell(length(contrast.matrix),size(y,2));
 results.GLM(3).F=zeros(length(contrast.matrix),size(y,2));
 results.GLM(3).Fsig=zeros(length(contrast.matrix),size(y,2));
 results.GLM(3).ces=zeros(length(contrast.matrix),size(y,2));
 results.GLM(3).ces_vector=cell(length(contrast.matrix),size(y,2));
 results.GLM(3).dof=zeros(length(contrast.matrix),size(y,2));           

  for i=1:size(y,2)
     disp(['GLM(3) ' num2str(i)])
     temp_y=y(:,i);
     if rank(X3(~isnan(temp_y),:)'*X3(~isnan(temp_y),:))<rank(X3)
         continue
     end

     [beta, rvar, ~, ~] = fast_glmfit(temp_y(~isnan(temp_y)),X3(~isnan(temp_y),:));
     
     for c=1:length(contrast.matrix)
         [F, Fsig, ces, edof] = fast_fratio(beta,X3(~isnan(temp_y),:),rvar,contrast.matrix{c});
         ces_covar=rvar*contrast.matrix{c}*crazy_mat3*contrast.matrix{c}';
         ces_std=sqrt(mean(diag(ces_covar)))/sqrt(length(diag(ces_covar)));
         results.GLM(3).ces_std(c,i)=ces_std;
         results.GLM(3).ces_covar{c,i}=ces_covar;
         results.GLM(3).F(c,i)=F;
         results.GLM(3).Fsig(c,i)=Fsig;
         results.GLM(3).ces(c,i)=mean(ces);
         results.GLM(3).ces_vector{c,i}=ces;
         results.GLM(3).dof(c,i)=edof;
     end
     
     
     results.GLM(3).beta(:,i)=beta;
     results.GLM(3).rvar(1,i)=rvar;
     
 end
 results.GLM(3).X=X3;
 results.GLM(3).contrast.name=contrast.name;
 results.GLM(3).contrast.matrix=contrast.matrix;

%To check ortogonality of the matrix:
% sx=X3.^2;
% a=(sum(sx,1)).^(.5);
% normat=repmat(a,size(X3,1),1);
% X0=X3./normat;
% imagesc(X0'*X0)