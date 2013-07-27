function v = stdv(X, xmean)

if nargin == 1
    xmean = mean(X);
end

diff = bsxfun(@minus, X, xmean);
sq_dist = sum(diff.^2, 2);
v = sqrt( mean( sq_dist ) );