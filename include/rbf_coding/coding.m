function code = coding(model, data, knn)

len = size(data, 1);
 
dist = dist2(data, model.dict); 
[~, help_ind] = sort(dist, 2);
 
ind = (help_ind(:, 1:knn)-1)*len + repmat([1:len]', 1, knn);
sigma2 = repmat( model.sigma', len, 1 );
sigma2 = sigma2(ind);

code = zeros(size(dist), 'single');
code(ind) = exp( -0.5*dist(ind)./(sigma2.^2) ) ./ ( sigma2*sqrt(2*pi) );    % calculate code using gaussian function

code = bsxfun(@rdivide, code, sum(code,2));     % L1 normalization



function n2 = dist2(x, c)
%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
  		ones(ndata, 1) * sum((c.^2)',1) - ...
  		2.*(x*(c'));
