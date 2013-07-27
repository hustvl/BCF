% The algorithm:
%	
%	1. Initialize the graph
%		for each p1
%			for each p2
%				if (p1,p2) inside the shape
%					add edge e(p1,p2) into E
%					set dis_mat(p1,p2)=dis_mat(p2,p1)
%					set ang_mat(p1,p2), ang_mat(p2,p1)
%					set viewable(p1,p2)=viewable(p2,p1)=1
%
%	2. for each p1, 
%		2.1 compute the distance and angle from it to all other points
%			initialize edge list E1
%			"bellman_ford"
%		2.2 compute shape context at p1

function [sc,V,E,dis_mat,ang_mat] = compu_contour_innerdist_SC( ...
											cont, fg_mask, ...
											n_dist, n_theta, ...
											bTangent, bSmoothCont, ifig)

if ~exist('bSmoothCont')	bSmoothCont = 0;	end
if ~exist('bTangent')		bTangent = 0;		end

%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);
V	= cont;


%- the following build the graph, in that each pair of points has an edge
%between them if they can see each other inside the shape boundary
E	= build_graph_contour_C(X,Y,fg_mask,bSmoothCont);
E	= E';

% disp_graph(V,E);		keyboard;

%-- the following compute the shortest path from the above graph, the
%length of the shortest path is used as the inner-distance. At the same
%time, the ang_mat shows the angles which is used as the "inner-angle" for
%building the inner-distance shape context
[dis_mat,ang_mat] = bellman_ford_ex_C(X,Y,E);


%-- Normalize shape context with the tangent orientation, this is to
%compute the inner-angle
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
n_pt		= size(V,1);
dists		= zeros(n_pt-1,n_pt);		% distance and orientation matrix, the i-th COLUMN contains
angles		= zeros(n_pt-1,n_pt);		% dist and angles from i-th point to all other points
id_gd		= setdiff(1:n_pt*n_pt, 1:(n_pt+1):n_pt*n_pt);
dists(:)	= dis_mat(id_gd);
angles(:)	= ang_mat(id_gd);

% tic;
[sc]		= comp_sc_hist(dists,angles,n_dist,n_theta);
% disp(['computing shape context take ' i2s(toc) ' seconds']);

%-- demonstrate
if exist('ifig') & ifig>0
	figure(ifig);	clf; hold on;	imagesc(fg_mask);colormap(gray);	axis equal;	
	sGraph	= ['|V|=' i2s(size(V,1)) ', |E|=' i2s(size(E,1))];
	disp_graph(V,E);		title(sGraph);
	
	figure(ifig+4);
	disp_schist(fg_mask,V,sc,n_dist,n_theta);
	
	keyboard;
end

return;
