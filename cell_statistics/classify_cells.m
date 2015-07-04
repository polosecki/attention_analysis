xdata=pca_input_matrix;
group=is_pitd;

max_perm=50;
%posta=0.3113;

c = cvpartition(length(group),'leaveout');
%c = cvpartition(length(group),'kfold',round(.9*length(group)));
posta=crossval('mcr',xdata,group,'Predfun',@svm_decoding,'partition',c);

for i=1:max_perm
disp(['running permutation number ' num2str(i)])
c = cvpartition(length(group),'leaveout');
mcr(i) = crossval('mcr',xdata,group(randperm(length(group))),'Predfun',@svm_decoding,'partition',c);%'holdout',0.05,'kfold',20)
end