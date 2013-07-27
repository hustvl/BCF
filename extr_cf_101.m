
function extr_cf_101(img_dir, fmt, sc_para)
    
    database = retr_database_dir(img_dir, fmt);

	tic
    % matlabpool(4);
    parfor n = 6908:length(database.path)
        fprintf('Extracting contour fragment: %d of %d\n', n, length(database.path));
        func_extr_cf(database.path{n}, sc_para ); 
    end
    % matlabpool close;
	toc

function func_extr_cf(im_name, cf_para )
    
	I = imread( im_name );
    I = im2bw( I, cf_para.t_edge );
    
    [yy, xx] = find( I > 0 );
    miny = max(min(yy)-1,1);
    minx = max(min(xx)-1,1);
    maxy = min(max(yy)+1,size(I,1));
    maxx = min(max(xx)+1,size(I,2));
    I = I( miny:maxy, minx:maxx );
    
    
    	I = bwmorph(I, 'thin', inf);
    	els = edgelink( I, 10 ); 		% disp(length(els));
	els = relinkedge_iter(els, 10, 5, 0.6);	% disp(length(els));
    
    pnts = {};
    for i = 1:length(els)
        C = els{i};
        C(:,1) = els{i}(:,2);
        C(:,2) = els{i}(:,1);
        C = C(1:end-1, :);

	if size(C,1) < cf_para.min_len_cf
		continue;
	end
		
 
        [XIs, YIs] = uniform_interp( C(:,1), C(:,2), round( size(C, 1) * 0.99 )-1);
        C = [C(1,:); [XIs YIs]];
        
        pnts_ = {};
        if size(C,1) > cf_para.n_pntsamp 
            pnts_ = extr_raw_pnts( C, cf_para.max_curvature, cf_para.n_contsamp, cf_para.n_pntsamp );
        else
            if size(C,1) > cf_para.min_len_cf
                pnts_ = cell(1,1);
                pnts_{1} = C;
            end
        end
        
        pnts = [pnts; pnts_];
    end
     
    
%     if strcmp( im_name( end-8:end-7 ), '10' )
%         figure('color', 'w'),
%         for i = 1:length(pnts)   
%             ps = pnts{i};
%             subplot( ceil(length(pnts)/10), 10, i);
%             plot(ps(:,1), ps(:,2));
%             axis off;
%             axis equal;
%         end
%         saveas( gcf, ['fig\', im_name( end-8:end-4), '.png'] );
%         close all;
%     end
        
    len = length(pnts);

    feat_sc = zeros( 0, len ); 
    xy = zeros( len, 2 );
    
    for i = 1:len
	cf = pnts{i};	

	if ~is_clockwise(cf)
		cf = cf(end:-1:1, :);
	end

        sc = shape_context( cf );
        sc = sc(:);
        sc = sc / sum(abs(sc));
        
        feat_sc(  1:size(sc,1), i ) = sc;

        xy( i, 1:2 ) = mean(cf); %sc( round(end/2), : );
    end
    
    sz = size(I);
    save([im_name(1:end-4), '_feat.mat' ], 'feat_sc', 'pnts', 'xy', 'sz');
    
    
