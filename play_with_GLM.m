[results]=make_GLM_fun(29,'Quincy','PITd','time_course',0,0);

%%
mat_used=1;
t_used=1; % in seconds
[~,t_idx]=min(abs(results{mat_used}.time-t_used));

y=results{mat_used}.y(:,t_idx);

GLM_used=3;
X=results{mat_used}.GLM(GLM_used).X;
linear_betas=results{mat_used}.GLM(GLM_used).beta(:,t_idx);
[b,dev,stats] = glmfit(X,y,'poisson','link','identity','constant','off','estdisp','off');
[b,dev,stats] = glmfit(X,y,'normal','link','identity','constant','off','estdisp','off');
linear_se=results{mat_used}.GLM(GLM_used).ces_std(:,t_idx)
[stats.se/stats.s linear_se]
[stats.p linear_p]

%%

        contrasts_plotted={logical([1 0 0 0]);
                           logical([0 1 0 0]);
                           logical([1 1 0 0 0 0 0 1 1 0 1 0 1 0 0 0])};
%                           logical([1 1 1 0 1 0 0 1 1 0 1 0 1 0 0 0])};                       
plot_GLM_contrasts(results,contrasts_plotted,0,'Quincy_PITd',29);
