function [res_dist,dis_ac,dis_sc,dis_be,match_cost,n_match] = dist_bw_sc( ...
											sc1,sc2, V1,V2, n_dist,n_theta, ...
											im1,im2, w, dis_sc, ifig)

if ~exist('dis_sc') 	dis_sc=[];	end
if ~exist('ifig') 		ifig=-1;	end

%-- compute SC distance b/w sc1 and sc2
% if isempty(dis_sc)		dis_sc	= dist_bw_sc_C( sc1, sc2, 0);	end
[dis_sc,costmat]	= dist_bw_sc_C( sc1, sc2, 0);

dis_be		= 0;
dis_ac		= 0;
match_cost	= 0;
n_match		= 0;
% % calculate shape context cost
% [a1,b1]=min(costmat,[],1);
% [a2,b2]=min(costmat,[],2);
% sc_cost=max(mean(a1),mean(a2));

%-- compute bipartite matching

%- Compute cost matrix
n1			= size(sc1,2);
n2			= size(sc2,2);
b_add_dummy	= 0;
if ~b_add_dummy
	costmat2	= costmat';
else
	% add dummy points
	eps_dum		= 0.25;
	dum_ratio	= 0;

	n_new		= round(max(n1,n2)*(1+dum_ratio));
	costmat2	= eps_dum*ones(n_new,n_new);
	costmat2(1:n2,1:n1)	= costmat';
end

%- MATCHING
% tic;
% bCircularMatching	= 0;
if 1
% 	tic;	[cvec1,match_cost1]	= DPMatching(costmat2,bCircularMatching);	toc;tic
	thre	= .5 * mean(mean(costmat2));
	search_step	= 1;
	num_start	= 2;
	[cvec,match_cost]	= DPMatching_C(costmat2,thre,num_start,search_step);
% 	toc,	keyboard
elseif 1
	tic;	[cvec,match_cost]	= hungarian(costmat2);		toc
else
%	[rowSol, colSol] = lap_mt( costmat2 );
	n_pt	= size(costmat2,1);
	sData	= 'costmat.dat';
	pF		= fopen(sData,'w');
	sFormat	= [repmat('%f  ', 1, n_pt), '\n'];
	fprintf(pF,sFormat,costmat2);
	fclose(pF);
	
	sRes	= 'cvec.txt';
	scmd	= ['lap_test ' sData ' ' sRes ' ' i2s(n_pt)];
	dos(scmd);
	cvec	= load(sRes);
	cvec	= cvec(2,:);
	
	if cvec(1)==0		[cvec,cst]	= hungarian(costmat2);		end	
% 	disp(toc)
end
% toc
% keyboard

% extract coordinates of non-dummy correspondences and use them to estimate transformation
id_gd1		= find(cvec<=n2);
id_gd1		= id_gd1(find(id_gd1<=n1));
id_gd2		= cvec(id_gd1);

pt_from		= V1(id_gd1,:);
pt_to		= V2(id_gd2,:);
n_match		= length(pt_from);

if w(3)>.00001 & 1
% 	keyboard
	beta_k		= 0;
	[cx,cy,E]	= bookstein(pt_from,pt_to, beta_k);
	dis_be		= E;
end

res_dist	= w(1)*dis_ac + w(2)*dis_sc + w(3)*(dis_be/n_match) + w(4)*match_cost;

% display correspondence
% ifig	= 12;
if exist('ifig') & ifig>0
	figure(ifig); clf;	hold on;
	plot(V1(:,1),V1(:,2),'r.');
	plot(V2(:,1)+150,V2(:,2)+150,'b.');
	
	plot([pt_from(:,1) pt_to(:,1)+150]', [pt_from(:,2) pt_to(:,2)+150]', '+-');

	title(['#match=' i2s(n_match) ' cost=' num2str(match_cost)]);
	xlabel(['dist=' num2str(res_dist) ' ac=' num2str(dis_ac) ' sc=' num2str(dis_sc) ' be=' num2str(dis_be)]);
	keyboard
end

return;

