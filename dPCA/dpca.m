% DISCLAIMER
% ---------------------------------------------------------------------
% This work is released under the
%
% Creative Commons Attribution-NonCommercial 3.0 Unported (CC BY-NC 3.0)
%
% license. Therefore, you are free to copy, redistribute and remix
% the code. You may not use this work for commercial purposes (please 
% contact the authors at wieland.brendel@neuro.fchampalimaud.org).  
% You are obliged to reference the work of the original authors:
%
% Wieland Brendel & Christian Machens, published at NIPS 2011 "Demixed
% Principal Component Analysis", code@http://sourceforge.net/projects/dpca/
%
% USAGE AT YOUR OWN RISK! The authors may not be hold responsible for any
% kind of damages or losses of any kind that might be traced back to the
% usage or compilation of this work.
% ---------------------------------------------------------------------
    
function W = dpca(Y,comps,maxstep,tolerance)
    % Performs a DPCA analysis on the data set Y. For further information
    % on this method please check the paper
    %
    % http://books.nips.cc/papers/files/nips24/NIPS2011_1440.pdf
    %
    % or the manual/code at
    % 
    % http://sourceforge.net/projects/dpca/
    %
    % INPUT
    % -------
    %  Y: multidimensional array with the first index being the
    %  observed object (e.g. neuron number) and subsequent dimensions
    %  referring to different parameters. E.g. to access neuron 5 at time
    %  t=10 and stimulus=2 you write
    % 
    %  Y[5,10,2]
    %
    %  comps:   Number of latent dimensions (i.e. # of components)
    %  maxstep: maximum number of steps (default 250)
    %  tolerance: minimum relative change of the objective function
    %
    % RETURNS
    % -------
    %  W: loading matrix
    
    % load default parameters
    if isempty(tolerance), tolerance = 1d-7; end
    if isempty(maxstep), maxstep = 100; end
    
    % load covariance matrices    
    [covs, C] = dpca_covs(Y);
    covs = values(covs);
    
    % init loading matrix with PCA solution
    [W,D] = eigs(C,comps);
       
    % set parameters for line search
    t=1/2; a = 2/t;  % initial discount, step size & loop iterator
    lam = 2;		    % trade-off between variance and demixing (higher values = higher demixing)
    old_L = L(W,covs,lam);
    fprintf('value of b is %1.7e\n',old_L)
    steps = 0;
    xchange = 1;
           
    while steps < maxstep && xchange > tolerance
        Q = W;    
        Gw = Lgradfull(W,covs,C,lam);
              
        % Line search using quadratic interpolation
        % Wolfe conditions are hard to implement because the objective
        % is costly to evaluate (almost as costly as the gradient)
        
        % sample the objective at for three different step-sizes
        x = [0,a,2*a];
        x = x';
        y = zeros(3,1);
        y(1) = old_L;
        
        for i=2:3
            T = Q + x(i)*Gw;
            T = qn(T);
            y(i) = L(T,covs,lam);
        end

        step = 0;
        while all((y - y(1))./y < tolerance)
            % fprintf('All sample points below current objective. Adding samples')
            step = step*t;            
            x(size(x,1)+1) = step;

            T = Q + step*Gw;
            T = qn(T);
            y(size(y,1)+1) = L(T,covs,lam);
            
            if size(y,1) > 20
                fprintf('Maximum search depth reached.')
                break
            end
        end
               
        % fit a quadratic function to the three points        
        p = polyfit_without_warning(x,y,2);
        
        % the minimum of the quadratic function is set as the new step-size
        quad_step = -p(2)/(2*p(1));
                   
        W = Q + quad_step*Gw;
        W = qn(W);
        new_L = L(W,covs,lam);
        xchange = (new_L - old_L)/new_L;

        a = quad_step/2;
        old_L = new_L;
        steps = steps + 1;
        
        disp(['iteration ', num2str(steps), ' @ objective ', num2str(new_L)])
    end
    
    Wb = bsxfun(@times, W, 1./sqrt(sum(W.^2, 1)));
    for ll=1:size(covs,2)
       diag(Wb.'*covs{ll}*Wb)./diag(Wb.'*C*Wb);
    end
    
    function l = L(W,covs,lam)
        % evaluate loss function
        l = 0;
        % normalize the rows of W
        Wb = bsxfun(@times, W, 1./sqrt(sum(W.^2, 1)));
        % multiply with marginal cov matrices
        Y = zeros(size(W,2),size(covs,2));
        for ll=1:size(covs,2)
           Y(:,ll) = diag(Wb.'*covs{ll}*Wb);
           %Y(:,ll) = sum(Wb.'*covs{ll}.*Wb.',2);
        end
        Y = sqrt(Y);

        % calculate objective
        for kk=1:size(W,2)
            v = Y(kk,:);
            l = l + norm(v)^2*(norm(v)/norm(v,1))^lam;
        end
    end

    function Gw = Lgradfull(W,covs,S,lam)
        % return gradient of loss function
	    x1 = zeros(size(W,2),1);
	    dX2 = S*W;
	    x2 = sqrt(diag(W.'*dX2));
	    dX1 = zeros(size(W));
	
        for j=1:length(covs)
	      CkW = covs{j}*W;
	      varC = sqrt(diag(W.'*CkW));
	      x1 = x1 + varC;
          dX1 = dX1 + bsxfun(@rdivide, CkW, varC');
        end
        
        v = (x2./x1).^lam;
        v2 = x2.^(2+lam)./(x1.^(lam+1));
        
        Gw = (2+lam)*bsxfun(@times, dX2, v') - lam*bsxfun(@times, dX1, v2');
    end

    function W = qn(W)
        % symmetric orthogonalization
        W = W/norm(W);
        N = size(W,2);
        T = W.'*W;
        while norm(T - eye(N)) > 0.000000001
           W = W*(3*eye(N) - T)/2;
           T = W.'*W;
        end       
    end  

end
