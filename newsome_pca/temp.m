align=2
nPCAs_used=12;
cc=coeff(:,1:nPCAs_used);
% glm_betas=reshape(betas_population{align},size(betas_population{align},1),size(betas_population{align},2)*size(betas_population{align},3))';
% denoised_betas=glm_betas*(cc*cc');
% figure;
% subplot(1,2,1)
% imagesc(glm_betas,[-.1 .9])
% subplot(1,2,2)
% imagesc(denoised_betas,[-.1 .9])

%original_glm=betas_population{align};
%denoised_glm=reshape(denoised_betas,size(betas_population{align}));

denoised_betas_population=betas_population;
figure;
beta_directions{align}=zeros(size(betas_population{align},1),size(betas_population{align},3));
for beta_index=1:size(betas_population{align},3)
    beta_temp=squeeze(betas_population{align}(:,:,beta_index));
    denoised_temp=(cc*cc')*beta_temp;%beta_temp*(cc*cc');
    denoised_betas_population{align}(:,:,beta_index)=denoised_temp;
    beta_norm=sqrt(sum(denoised_temp.^2,1));
    [~,max_ind]=max(beta_norm);
    beta_directions{align}(:,beta_index)=denoised_temp(:,max_ind)/beta_norm(max_ind);
    cos_angle=beta_directions{align}(:,beta_index)'*(denoised_temp./repmat(beta_norm,37,1));
    subplot(2,2,beta_index)
    plot(100*sin(acos(cos_angle/max(cos_angle))))
    title(num2str(beta_index))
end
% Orthogonalize the vector base using QR decomposition:


show_index=3;
figure;
subplot(1,2,1)
imagesc(squeeze(betas_population{align}(:,:,show_index)),[-.1 .9])
subplot(1,2,2)
imagesc(squeeze(denoised_betas_population{align}(:,:,show_index)),[-.1 .9])

