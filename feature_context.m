function V = feature_context(features, fc_para)
%================================================
% 
% Usage:
% Combining Sparse Coding and Feature Context together. 
%
% Inputss:
% features        -structure defining the feature set of an image   
%                   .code     	encoded image feature
%                   .x          x locations of each local feature, 2nd
%                               dimension of the matrix
%                   .y          y locations of each local feature, 1st
%                               dimension of the matrix
%                   .wid	    width of the image
%                   .hgt	    height of the image
% fc_para		-strcture defining the parameter of Feature Context.
% 					.vertex     vertexes of the polar coordinates, reference to the image size
% 					.bin_num	a 1*2 matrix, [number_of_dists, number_of_orientation]
% 					.r_inner	
% 					.r_outter
% Output:
% 

vertex = fc_para.vertex;
bins_num = fc_para.bins_num;
r_inner = fc_para.r_inner;
r_outer = fc_para.r_outer;

delta_r = 0.05;
delta_theta = 0.05;

V = [];
for iter2 = 1:size(vertex, 2)

	dd = [ features.x-features.wid*vertex(1,iter2), features.y-features.hgt*vertex(2,iter2) ];
	
	% distance
	r = sqrt( dd(:, 1).^2 + dd(:, 2).^2 );
	r = r ./ mean(r);
	bin_r = logspace( log10(r_inner), log10(r_outer), bins_num(1)-1 );
    bin_r = [0, bin_r, inf];
	
	r_label  = zeros(1, length(r) );
	for n = 1 : length(bin_r)-1
		help_ind =  r>=(bin_r(n)-delta_r) & r<(bin_r(n+1)+delta_r) ;
		r_label(help_ind) = n;
    end
	
	% angle
	theta = atan2( dd(:, 1), dd(:, 2) ) + pi;
	theta = rem(rem(theta, 2*pi) + 2*pi, 2*pi);
	bin_theta = linspace(0,2*pi, bins_num(2)+1 );
	theta_label = zeros( 1, length(r) );
	for n = 1 :length(bin_theta)-1
		help_ind =  theta>=(bin_theta(n)-delta_theta) & theta<=(bin_theta(n+1)+delta_theta) ;
		theta_label(help_ind) = n;
	end

	space_label = sub2ind( bins_num(1:2), r_label, theta_label );
	
	grid_v = zeros( size(features.code,1), bins_num(1)*bins_num(2) );
	for iter3 = 1:bins_num(1)*bins_num(2)
        help_ind = find(space_label == iter3);
        if isempty(help_ind),
            continue;
        end   
        % max-pooling
		grid_v(:, iter3) = max( features.code(:, help_ind ), [], 2);
	end
	
	V = [V, grid_v(:)'];
end

 V = [V, max(features.code, [], 2)'];









