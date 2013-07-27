
function normalize_shape(img_dir, fmt)
    
    database = retr_database_dir(img_dir, fmt);

	tic
    % matlabpool(4);
    for n = 1:length(database.path)
        fprintf('Normalize shape: %d of %d\n', n, length(database.path));
        func_normalize_shape(database.path{n} ); 
    end
    % matlabpool close;
	toc

function func_normalize_shape(im_name )
    
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
    
    
        
    I = no_border_img(I, 1);
    I = imresize( I, sqrt(200000/(size(I,1)*size(I,2))) );
    I = no_border_img(I, 5);
    
    
    imwrite(I, [im_name(1:end-4), '.png']);
    
    
    
       
function bw2 = no_border_img(bw,s)
    [oy, ox] = find(bw>0);
    maxy = max(oy);
    miny = min(oy);
    maxx = max(ox);
    minx = min(ox);
    bw2 = zeros( maxy-miny+2*s, maxx-minx+2*s );
    bw2( s+1:maxy-miny+s+1, s+1:maxx-minx+s+1 ) = bw( miny:maxy, minx:maxx );
    bw2 = bw2 > 0;
    
    
    
    