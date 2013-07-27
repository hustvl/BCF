% ========================================================================
% USAGE: [Coeff]=LLC_coding_appr(B,X,knn,lambda)
% Approximated Locality-constraint Linear Coding
%
% Inputs
%       B       -M x d codebook, M entries in a d-dim space
%       X       -N x d matrix, N data points in a d-dim space
%       knn     -number of nearest neighboring
%       lambda  -regulerization to improve condition
%
% Outputs
%       Coeff   -N x M matrix, each row is a code for corresponding X
%
% Jinjun Wang, march 19, 2010
% ========================================================================

function [Coeff] = soft_coding(B, X, knn)

if ~exist('knn', 'var') || isempty(knn),
    knn = 5;
end


nframe=size(X,1);
nbase=size(B,1);

% find k nearest neighbors
XX = sum(X.*X, 2);
BB = sum(B.*B, 2);
D  = repmat(XX, 1, nbase)-2*X*B'+repmat(BB', nframe, 1);

% IDX = zeros(nframe, knn);
% for i = 1:nframe,
% 	d = D(i,:);
% 	[dummy, idx] = sort(d, 'ascend');
% 	IDX(i, :) = idx(1:knn);
% end

[Ds, IDX] = sort(D, 2); 
IDX = IDX(:,1:knn);

% llc approximation coding
Coeff = zeros(nframe, nbase);
for i=1:nframe
   idx = IDX(i,:);
   dis = D(i,idx);
   
   w = exp( -dis/mean(dis) );
   w = w  / sum(w);
   
   Coeff(i,idx) = w;
end
