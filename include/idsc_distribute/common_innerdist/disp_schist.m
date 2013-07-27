function disp_schist(im,pts,sc_hist,n_dist,n_theta,nRand)
	
	if ~exist('nRand')		nRand = -20;	end
	
	N	= size(pts,1);
	X	= pts(:,1);
	Y	= pts(:,2);
	
	if nRand>0
		for ii=1:round(nRand)
			v		= ceil(round(rand*N));
			sctmp	= reshape(sc_hist(:,v),n_dist,n_theta);
			subplot(1,2,1);	hold on;
			imagesc(sctmp); colorbar;	%colormap(gray); 
			title(i2s(v));	xlabel('[-pi,pi]');	drawnow
			
			subplot(1,2,2);	hold on;
			imshow(im);	axis equal; axis xy; axis on;
			plot(X,Y,'b.');
			plot(X(v),Y(v),'r*');	drawnow;
			pause;
		end
	else
		for v=1:2:N
			sctmp	= reshape(sc_hist(:,v),n_dist,n_theta);
			subplot(1,2,1);	hold on;
			imagesc(sctmp); colorbar;	%colormap(gray); 
			title(i2s(v));	xlabel('[-pi,pi]');	drawnow
			
			subplot(1,2,2);	hold on;
			imshow(im);	axis equal;	axis xy; axis on;
			plot(X,Y,'b.');
			plot(X(v),Y(v),'r*');	drawnow;
			pause;
		end
	end
	
return;
