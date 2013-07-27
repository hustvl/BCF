% th
function [sc,V,dis_mat,ang_mat] = shape_context( cont, n_ref, n_dist, n_theta, bTangent)

if ~exist('n_dist')		n_ref = 5;		end			
if ~exist('n_dist')		n_dist = 5;		end											
if ~exist('n_theta')		n_theta = 12;		end											
if ~exist('bTangent')		bTangent = 1;		end

%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);

%------ Set reference point
si = round( linspace(1, n_pt, n_ref) );
V	= cont( si, : );
vx  = V(:, 1);
vy  = V(:, 2);

%-- Orientations and geodesic distances between all landmark points ---------
Xs			= repmat(X,1,n_ref);
vxs			= repmat(vx',n_pt,1);
dXs			= Xs-vxs;
Ys			= repmat(Y,1,n_ref);
vys			= repmat(vy',n_pt,1);
dYs			= Ys-vys;
dis_mat		= sqrt(dXs.^2+dYs.^2);
ang_mat		= atan2(dYs,dXs);


%-- Normalize shape context with the tangent orientation --------------------
if bTangent
	Xs	= [X(end); X; X(1)];
	Ys	= [Y(end); Y; Y(1)];

	gX	= gradient(Xs);
	gY	= gradient(Ys);
	thetas	= atan2(gY,gX);
	thetas	= thetas(2:end-1);
	thetas	= repmat(thetas,1,n_ref);
	
	ang_mat	= ang_mat-thetas;
	idx		= find(ang_mat>pi);
	ang_mat(idx)	= ang_mat(idx)-2*pi;
	idx		= find(ang_mat<-pi);
	ang_mat(idx)	= ang_mat(idx)+2*pi;
end


%-- Compute shape context  ----------------------------------------
% dists		= zeros(n_pt-1,n_pt);		% distance and orientation matrix, the i-th COLUMN contains
% angles		= zeros(n_pt-1,n_pt);		% dist and angles from i-th point to all other points
% id_gd		= setdiff(1:n_pt*n_pt, 1:(n_pt+1):n_pt*n_pt);
% dists(:)	= dis_mat(id_gd);
% angles(:)	= ang_mat(id_gd);

[sc]		= shape_context_core(dis_mat,ang_mat,n_dist,n_theta);

return;





function [sc_hist] = shape_context_core(dists,angles,n_dist,n_theta)

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
	
	sc_hist	= sc_hist/size(dists,1);
	
return;
