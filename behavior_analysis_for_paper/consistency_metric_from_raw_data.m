function [success, chance, confidence, p_val]=consistency_metric_from_raw_data(RT_absolute,BRT_time,nBRT_time,sacc_dir,BRT_dir,nBRT_dir,n_perms)



detect_comp=(RT_absolute-BRT_time)>0 & (RT_absolute-BRT_time)<2.7;
detect_ncomp=(RT_absolute-nBRT_time)>0 & (RT_absolute-nBRT_time)<2.7;

disc_comp = nan(size(detect_comp));
disc_comp(detect_comp & (sacc_dir==BRT_dir))=true;
disc_comp(detect_comp & (sacc_dir~=BRT_dir))=false;
disc_ncomp = nan(size(detect_ncomp));
disc_ncomp(detect_ncomp & (sacc_dir==nBRT_dir))=true;
disc_ncomp(detect_ncomp & (sacc_dir~=nBRT_dir))=false;


%detect.both=mean(detect_comp & detect_ncomp)*100;
%detect.targ=mean(detect_comp & ~detect_ncomp)*100;
%detect.dist=mean(~detect_comp & detect_ncomp)*100;
%detect.none=mean(~detect_comp & ~detect_ncomp)*100;

success.both=mean(detect_comp & detect_ncomp & disc_comp==1 & disc_ncomp==1)*100;
success.targ=mean(detect_comp & disc_comp==1 & (~detect_ncomp | disc_ncomp~=1))*100;
success.dist=mean(detect_ncomp & disc_ncomp==1 & (~detect_comp | disc_comp~=1))*100;
success.none=mean((~detect_ncomp | disc_ncomp~=1) & (~detect_comp | disc_comp~=1))*100;

dirs=unique(BRT_dir);
t{1}=find(ismember(BRT_dir,[dirs(1) dirs(3)]));
t{2}=find(ismember(BRT_dir,[dirs(2) dirs(4)]));

normal_idx=[t{1};t{2}];

distr.both=zeros(n_perms,1);
distr.targ=zeros(n_perms,1);
distr.dist=zeros(n_perms,1);
distr.none=zeros(n_perms,1);

for i=1:n_perms
    rand_idx=[t{1}(randperm(length(t{1}))); t{2}(randperm(length(t{2})))];
    
    detect_comp=(RT_absolute(normal_idx)-BRT_time(rand_idx))>0 & (RT_absolute(normal_idx)-BRT_time(rand_idx))<2.7;
    detect_ncomp=(RT_absolute(normal_idx)-nBRT_time(rand_idx))>0 & (RT_absolute(normal_idx)-nBRT_time(rand_idx))<2.7;
    
    disc_comp = nan(size(detect_comp));
    disc_comp(detect_comp & (sacc_dir(normal_idx)==BRT_dir(rand_idx)))=true;
    disc_comp(detect_comp & (sacc_dir(normal_idx)~=BRT_dir(rand_idx)))=false;
    disc_ncomp = nan(size(detect_ncomp));
    disc_ncomp(detect_ncomp & (sacc_dir(normal_idx)==nBRT_dir(rand_idx)))=true;
    disc_ncomp(detect_ncomp & (sacc_dir(normal_idx)~=nBRT_dir(rand_idx)))=false;
    distr.both(i)=mean(detect_comp & detect_ncomp & disc_comp==1 & disc_ncomp==1)*100;
    distr.targ(i)=mean(detect_comp & disc_comp==1 & (~detect_ncomp | disc_ncomp~=1))*100;
    distr.dist(i)=mean(detect_ncomp & disc_ncomp==1 & (~detect_comp | disc_comp~=1))*100;
    distr.none(i)=mean((~detect_ncomp | disc_ncomp~=1) & (~detect_comp | disc_comp~=1))*100;
end

p_val.both = sum(abs(distr.both-mean(distr.both))>abs(success.both-mean(distr.both))) / n_perms;
p_val.targ = sum(abs(distr.targ-mean(distr.targ))>abs(success.targ-mean(distr.targ))) / n_perms;
p_val.dist = sum(abs(distr.dist-mean(distr.dist))>abs(success.dist-mean(distr.dist))) / n_perms;
p_val.none = sum(abs(distr.none-mean(distr.none))>abs(success.none-mean(distr.none))) / n_perms;

chance.both=median(distr.both);
chance.targ=median(distr.targ);
chance.dist=median(distr.dist);
chance.none=median(distr.none);

confidence.both = [prctile(distr.both,5) prctile(distr.both,95)];
confidence.targ = [prctile(distr.targ,5) prctile(distr.targ,95)];
confidence.dist = [prctile(distr.dist,5) prctile(distr.dist,95)];
confidence.none = [prctile(distr.none,5) prctile(distr.none,95)];

%todo: compute the distribution and take percentiles/p-values;
