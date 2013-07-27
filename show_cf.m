
function show_cf(img_dir, fmt, sc_para)
    
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
    
    I = double( I );
    C = extract_longest_cont(I, para.n_shapesamp);
    figure, imshow(I);
    hold on;
    plot( C(:,1), C(:,2), '.r', 'markersize', 5 );
    saveas(gcf, ['fig\cntr\', im_name(end-8:end-4), '_cntr.jpg' ] );
    close all
    
