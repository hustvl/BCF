
function extr_cf(img_dir, fmt, sc_para)
    
    database = retr_database_dir(img_dir, fmt);

	tic
    % matlabpool(4);
    parfor n = 1 : length(database.path)
        fprintf('Extracting contour fragment: %d of %d\n', n, length(database.path));
        func_extr_cf(database.path{n}, sc_para ); 
    end
    % matlabpool close;
	toc

function func_extr_cf(im_name, cf_para )
    
	I = imread( im_name );
    if size(I, 3) > 1
        for i=1:size(I,3)
            a=(I(:,:,i));
            if max(max(a))~=min(min(a))
                I=a;
                break               
            end
        end
    end
    
    I = double( I );
     
    if 0
        load([im_name(1:end-4), '_feat.mat' ], 'pnts');
    else
        C = extract_longest_cont(I, cf_para.n_shapesamp);
        
        sigma = 0.2;
        C = C + sigma * randn(size(C));
        
        pnts = extr_raw_pnts( C, cf_para.max_curvature, cf_para.n_contsamp, cf_para.n_pntsamp );
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
        sc = shape_context( pnts{i} );
        sc = sc(:);
        sc = sc / sum(sc);
        
        feat_sc(  1:size(sc,1), i ) = sc;
        xy( i, 1:2 ) = sc( round(end/2), : );
    end
    
    sz = size(I);
    save([im_name(1:end-4), '_feat.mat' ], 'pnts', 'feat_sc', 'xy', 'sz');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
