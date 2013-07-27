% [sc_hist]	= comp_sc_hist(dists,angles);
%		
% 	Compute shape context from input distances and angles
% Inputs: 	
% 		dists	: distance and orientation matrix, the i-th COLUMN contains
% 		angles	: dist and angles from i-th point to all other points

function [sc_hist] = comp_sc_hist(dists,angles,n_dist,n_theta)

	if ~exist('n_dist')		n_dist		= 10;	end
	if ~exist('n_theta')	n_theta		= 16;	end
	n_pts		= size(dists,2);

	%- preprocessing distances
	if 1
		% using log distance
		logbase		= 1.5;
		mean_dis	= mean(dists(:));
		b			= 1;
		a			= (logbase^(.75*n_dist)-b)/mean_dis;
		
		dists		= floor(log(a*dists+b)/log(logbase));
		dists		= max(dists,1);
		dists		= min(dists,n_dist);
	else
		% using linear distance
		logbase		= 1.5;
		mean_dis	= mean(dists(:));
		delta_dis	= mean_dis/n_dist;
		dists		= ceil(dists/delta_dis);
		dists		= max(dists,1);
		dists		= min(dists,n_dist);
	end		
% 	keyboard;
	
	%- preprocessing angles
	delta_ang	= 2*pi/n_theta;
	if 0
		angles		= angles+delta_ang/2;
		idx			= find(angles<0);
		angles(idx)	= angles(idx)+2*pi;
		angles		= ceil(angles/delta_ang);
	else
		angles		= ceil((angles+pi)/delta_ang);
	end
	angles		= max(angles,1);
	angles		= min(angles,n_theta);
	
	
	%- shape context
	sc_hist		= zeros(n_theta*n_dist, n_pts);
	sctmp		= zeros(n_dist,n_theta);
	for v=1:n_pts
		for dis=1:n_dist
			for ang=1:n_theta
				sctmp(dis,ang)	= length(find(dists(:,v)==dis & angles(:,v)==ang));
			end
		end
		sc_hist(:,v)	= sctmp(:);
		
		if 0
			figure(324);	clf; hold on;
			imagesc(sctmp); colorbar;	%colormap(gray); 
			title(i2s(v));	drawnow
			pause;
		end
	end
	
	sc_hist	= sc_hist/(n_pts-1);
	
return;