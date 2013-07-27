
function show_cf_occlusion(img_dir, fmt, sc_para)
    
    database = retr_database_dir(img_dir, fmt);

	tic
    for n = 1:length(database.path)
        fprintf('Show contour fragment: %d of %d\n', n, length(database.path));
        func_show_cf(database.path{n}, sc_para ); 
    end


function func_show_cf(im_name, para )
    
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
    
    figure, imshow(I);
    hold on;
    
    I = double( I );
    E = edge( I, 'sobel' );
    E = bwmorph(E, 'thin', inf);
    
    [ey, ex] = find( E > 0 );
    minx = min(ex);
    maxx = max(ex);
    t = round(minx + para.occlusion * (maxx - minx))-1;
    E( :, 1:t ) = 0;
    
    els = edgelink( E, 5 );
    
    colors = { '.r', '.b', '.y', '.k' };
    
    pnts = {};
    for i = 1:length(els)
        C = els{i};
        C(:,1) = els{i}(:,2);
        C(:,2) = els{i}(:,1);
        C = C(1:end-1, :);
        
        plot( C(:,1), C(:,2), colors{i}, 'markersize', 5 );
        
%         pnts_ = {};
%         if size(C,1) > para.n_pntsamp 
%             [XIs, YIs] = uniform_interp( C(:,1), C(:,2), para.n_shapesamp-1);
%             C = [C(1,:); [XIs YIs]];
%             pnts_ = extr_raw_pnts( C, para.max_curvature, para.n_contsamp, para.n_pntsamp );
%         else
%             if size(C,1) > 30
%                 pnts_ = cell(1,1);
%                 pnts_{1} = C;
%             end
%         end
%         
%         pnts = [pnts; pnts_];
    end

    
    saveas(gcf, ['fig\ed\', im_name(end-8:end-4), '_edge.jpg' ] );
    close all
    
