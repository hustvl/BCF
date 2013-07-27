addpath common_innerdist;
clear;

%------ Parameters ----------------------------------------------
ifig		= 1;
sIms		= {'01_01_m.bmp', '01_04_m.bmp'};

%-- shape context parameters
n_dist		= 5;
n_theta		= 12;
bTangent	= 1;
bSmoothCont	= 1;
n_contsamp	= 100;


%-- Extract inner-distance shape context
figure(ifig);	clf; hold on;	set(ifig,'color','w'); colormap(gray);
for k=1:2
	
	%- Contour extraction
	ims{k}	= double(imread(sIms{k}));
	Cs{k}	= extract_longest_cont(ims{k}, n_contsamp);
	
	%- inner-dist shape context
	msk		= ims{k};%>.5;
	[sc,V,E,dis_mat,ang_mat] = compu_contour_innerdist_SC( ...
									Cs{k},msk, ...
									n_dist, n_theta, bTangent, bSmoothCont,...
									0);
	scs{k}	= sc;

	%- demo code
	subplot(2,3,1+3*(k-1));	
	imagesc(ims{k});	hold on;
	sGraph	= ['graph for sh-path (zoom in for detail)'];%, |V|=' i2s(size(V,1)) ', |E|=' i2s(size(E,1))];
	disp_graph(V,E);		title(sGraph);

	subplot(2,3,2+3*(k-1));	
	ptids	= [1 25 50 75];
	plot(V(:,1),V(:,2),'b-');  axis equal; hold on;
	plot(V(ptids,1),V(ptids,2), 'xr');  
	axis off;
	title(['Four marked points and their IDSC (right)']);
	

	gid		= [5 6 11 12];
	sid		= {'bottom', 'right', 'top', 'left'};
	for p=1:length(ptids)
		v		= ptids(p);
		sctmp	= reshape(scs{k}(:,v),n_dist,n_theta);
		subplot(4,6,gid(p)+12*(k-1));
		imagesc(sctmp);
		title(['IDSC at pt ' i2s(v) ', ' sid{p}]);	
		xlabel('[-pi,pi]');	drawnow
		axis off;
	end
	
end


%-- compute SC distance b/w sc1 and sc2 and cost matrix
[dis_sc,costmat]	= dist_bw_sc_C( scs{1},scs{2}, 0);

%-- MATCHING
n_match		= 1;
num_start			= round(n_contsamp);
search_step			= 1;
thre				= .6;
[cvec,match_cost]	= DPMatching_C(costmat,thre,num_start,search_step);


%-- point correspondences for computation of 
n_pt		= size(scs{1},2);
id_gd1		= find(cvec<=n_pt);
id_gd2		= cvec(id_gd1);
pt_from		= Cs{1}(id_gd1,:);
pt_to		= Cs{2}(id_gd2,:);
n_match		= length(pt_from);

%-- display correspondence
ifig	= ifig+1;
figure(ifig); clf;	hold on; set(ifig,'color','w');
d	= 1;
plot([pt_from(1:d:end,1) pt_to(1:d:end,1)+150]', ...
	 [pt_from(1:d:end,2) pt_to(1:d:end,2)]', '+--',...
	 'linewidth',.5);

plot([Cs{1}(:,1); Cs{1}(1,1)],[Cs{1}(:,2); Cs{1}(1,2)],'k-', ...
	 [Cs{2}(:,1); Cs{2}(1,1)]+150, [Cs{2}(:,2); Cs{2}(1,2)],'k-',...
	 'linewidth',.5);
 
%plot(Cs{1}(1,1),Cs{1}(1,2),'rd', Cs{2}(1,1)+150,Cs{2}(1,2)-50,'bd');
sTitle	= ['#match=' i2s(n_match) ' cost=' num2str(match_cost)];
title(sTitle);
axis equal;	axis ij; axis off;

return;