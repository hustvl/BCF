% function disp_graph(V,E)

function disp_graph(V,E)
	
	N	= size(V,1);
	X	= V(:,1);
	Y	= V(:,2);
	Xs	= [X(E(:,1)),X(E(:,2))]';
	Ys	= [Y(E(:,1)),Y(E(:,2))]';
	plot(Xs,Ys,'r');
	plot(X,Y,'b.');

return;
