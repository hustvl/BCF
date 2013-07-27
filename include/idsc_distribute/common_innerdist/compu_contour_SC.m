% th
function [sc,V,dis_mat,ang_mat] = compu_contour_SC( cont, n_dist, n_theta, bTangent)
											
if ~exist('bTangent')		bTangent = 0;		end

%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);
V	= cont;

%-- Orientations and geodesic distances between all landmark points ---------
Xs			= repmat(X,1,n_pt);
dXs			= Xs-Xs';
Ys			= repmat(Y,1,n_pt);
dYs			= Ys-Ys';
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
	thetas	= repmat(thetas',n_pt,1);
	
	ang_mat	= ang_mat-thetas;
	idx		= find(ang_mat>pi);
	ang_mat(idx)	= ang_mat(idx)-2*pi;
	idx		= find(ang_mat<-pi);
	ang_mat(idx)	= ang_mat(idx)+2*pi;
end


%-- Compute shape context  ----------------------------------------
dists		= zeros(n_pt-1,n_pt);		% distance and orientation matrix, the i-th COLUMN contains
angles		= zeros(n_pt-1,n_pt);		% dist and angles from i-th point to all other points
id_gd		= setdiff(1:n_pt*n_pt, 1:(n_pt+1):n_pt*n_pt);
dists(:)	= dis_mat(id_gd);
angles(:)	= ang_mat(id_gd);

[sc]		= comp_sc_hist(dists,angles,n_dist,n_theta);

return;
