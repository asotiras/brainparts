function [W H] = opnmf_mem(X, K, w0, initMeth, max_iter, tol, ...
    iter0, save_step, outputdir)
%
% compute OPNMF based on Euclidean distance 
% Orthogonal Projective NonNegative Matrix Factorization
%
% The update rule is slightly modified to account for the high
% dimensionality of the data
%
% input:
%   X          nonnegative data input (D times N)
%   K          number of components
%   w0         given initialization of W (optional)
%   initMeth   determines which initialization method will be used if no w0
%              is provided
%   max_iter   maximum number of iterations (default 50000)
%   tol        convergence tolerance (default 1e-5)
%   iter0      initial iteration (used when resuming optimization after
%              possible failure - use in combination with saved 
%              intermediate results)
%   save_step  save intermediate results every # number of steps
%   outputdir  directory where intermediate results are saved
%
% output:
%   W          the factorizing matrix (D times K)
%   H          expansion coefficients   
%
% Relevant references:
%
% Sotiras, A., Resnick, S. M., & Davatzikos, C. (2015). Finding imaging 
% patterns of structural covariance via Non-Negative Matrix Factorization. 
% NeuroImage, 108, 1-16.
%
% Yang, Z., & Oja, E. (2010). Linear and nonlinear projective nonnegative 
% matrix factorization. Neural Networks, IEEE Transactions on, 21(5), 
% 734-749.

[Dinit,N] = size(X);

% Basic argument check
if ~isscalar(K) || ~isnumeric(K) || K<1 || K>min(Dinit,N) || K~=round(K)
    error('opnmf:badK','K should be positive integer no larger than the number of rows or columns in X');
end
if ~ismatrix(X) || ~isnumeric(X) || ~isreal(X) || any(any(X<0)) || any(any(~isfinite(X)))
    error('opnmf:badX','opnmf:X must be a matrix of non-negative values.')
end

% get rid of the background in order to boost the computational efficiency
% assumption: background corresponds to positions with stackwise mean value
% equal to zero
mean_im = mean(X,2);
data_matrix_nz = X((mean_im>0),:) ; 
X = data_matrix_nz ; clear data_matrix_nz ;

% variables
check_step = 100;
D = size(X,1);

% initialize w0
if ~exist('w0','var') || isempty(w0)
    disp('Initializing w0: ');
    if ~exist('initMet', 'var') || isempty(initMeth)
        initMeth = 1;
    end
    switch initMeth
        case 0
            disp('random initialization ...') ;
            W = rand(D,K);
        case 1 % NNDSVD
            disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
        case 2 % NNDSVDa
            disp('NNDSVDa initialization ...') ;
            [W,~] = NNDSVD(X,K,1) ;
        case 3 % NNDSVDar
            disp('NNDSVDar initialization ...') ;
            [W,~] = NNDSVD(X,K,2) ;
        case 4 % NNDSVD using randomized SVD
            disp('NNDSVD initialization using random SVD calculation for efficiency ...') ;
            [W,~] = NNDSVD(X,K,3) ;
        otherwise
            disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
    end
    disp('done') ;
else
    W = w0;
    clear w0 ;
end

% check variables
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = 50000;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
if (~exist('iter0','var') || isempty(iter0))
    iter0 = 1;
end
if (~exist('save_step','var') || isempty(save_step) )
    % save intermediate results (every 10% of the total process - default)
    save_step = floor(max_iter/10) ;
end
if (exist('outputdir','var') && ~isempty(outputdir))    
    % then create the directory
    if(~strcmp(outputdir(end),'/'))
        outputdir=[outputdir '/'];
    end
    if(~exist(outputdir,'dir'))
        success = mkdir(outputdir);
        if(~success)
            error('opnmf:BadDir',['Output directory ' outputdir ' can not be created']);
        end
    end
end

% start optimization
for iter=iter0:max_iter
    W_old = W;
    if mod(iter,check_step)==0
        fprintf('iter=% 5d ', iter);
    end
    
    % multiplicative update rule
    % the update rule is slightly modified to account for the high
    % dimensionality of the imaging data
    W = W .* (X*(X'*W)) ./ (W*((W'*X)*(X'*W)));    
    % As the iterations were progressing, computational time per iteration was increasing due to operations involving really small values
    W(W<1e-16)=1e-16;
    W = W ./ norm(W);
    
    diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    if diffW<tol
        fprintf('converged after %d steps.\n', iter);
        break;
    end
    
    if mod(iter,check_step)==0
        fprintf('diff=%.10f, ', diffW);
        fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
        fprintf('\n');
    end
    
    % save intermediate results
    if ( exist('outputdir','var') && ~isempty(outputdir) )
        if ( mod(iter,save_step) == 0 )
            fprintf('Saving intermediate results ...');
            save([outputdir 'IntermResultsExtractBases.mat'],'iter','D','K','W','-v7.3') ;
            fprintf('done\n') ;
        end
    end
end

% Reordering the output - putting them in standard form, loosely following
% what is done in the nnmf matlab function
H = W'*X ;

hlen = sqrt(sum(H.^2,2));
if any(hlen==0)
    warning(message('opnmf:LowRank', K - sum( hlen==0 ), K));
    hlen(hlen==0) = 1;
end
clear H
Wh = bsxfun(@times,W,hlen');

% Then order by w
[~,idx] = sort(sum(Wh.^2,1),'descend'); clear Wh
W = W(:,idx);
H = W'*X ;

% put results to original dimension
WW = zeros(Dinit,K) ;
WW(mean_im>0,:) = W ; clear W;
W = WW ; clear WW

