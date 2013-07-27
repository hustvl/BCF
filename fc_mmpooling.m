
function fc_mmpooling(img_dir, fea_dir)

database = retr_database_dir(img_dir, '*_code.mat');

fc_para.bins_num = [3, 6];
fc_para.r_inner = 0.2;
fc_para.r_outer = 1;

for  n = 1:length( database.path )
    label = database.label(n);
    path = database.path{n};
    load(path);
        
    fprintf('pyramid pooling %d of %d\n', n , length(database.path));
    
    im.code = code;
    im.wid = sz(1);
    im.hgt = sz(2);
    im.x = xy(:,1);
    im.y = xy(:,2);
    
    t = 3;
    [vx, vy] = meshgrid( [1/t:1/t:(t-1)/t]*sz(1), [1/t:1/t:(t-1)/t]*sz(2) );
    fc_para.vertex = [vx(:), vy(:)]';
    
    fea = feature_context(im, fc_para);
    fea = fea / sqrt(sum(fea.^2));
    
    [~, fname] = fileparts(path);
    feapath = fullfile(fea_dir, num2str(label), [fname(1:end-5) '.mat']);
    if ~isdir(fullfile(fea_dir, num2str(label))),
        mkdir(fullfile(fea_dir, num2str(label)));
    end  
    
    save(feapath, 'fea', 'label');
end


function feas = pyramid_pooling(pyramid, sz, xy, code)
    feas = zeros( 0, sum(pyramid.^2) );
    counter = 0;
    hgt = sz(1);
    wid = sz(2);
    x = xy(:,1);
    y = xy(:,2);
    code = code';
    
    for p = 1:length(pyramid) 
        for i = 1:pyramid(p)
            for j = 1:pyramid(p)
                
                yrang = hgt*[i-1,i]/pyramid(p);
                xrang = wid*[j-1,j]/pyramid(p);
                
                c_ = code( x >= xrang(1) & x <= xrang(2) & y >= yrang(1) & y <= yrang(2), : );
                f_ = max(c_);
                
                counter = counter + 1;
                feas( 1:length(f_), counter ) = f_;
            end
        end
    end
    
function feas = pascal_pooling(pyramid, sz, xy, code)
    feas = zeros( 0, 8 );
    counter = 0;
    hgt = sz(1);
    wid = sz(2);
    x = xy(:,1);
    y = xy(:,2);
    code = code';
    
    pyramid = [1, 2];
    
    for p = 1:length(pyramid) 
        for i = 1:pyramid(p)
            for j = 1:pyramid(p)
                
                yrang = hgt*[i-1,i]/pyramid(p);
                xrang = wid*[j-1,j]/pyramid(p);
                
                c_ = code( x >= xrang(1) & x <= xrang(2) & y >= yrang(1) & y <= yrang(2), : );
                f_ = max(c_);
                
                counter = counter + 1;
                feas( 1:length(f_), counter ) = f_;
            end
        end
    end

    yrang = hgt*[0,1];
    
    xrang = wid*[0,1/3];
    c_ = code( x >= xrang(1) & x <= xrang(2) & y >= yrang(1) & y <= yrang(2), : );
    f_ = max(c_);
    feas( 1:length(f_), 6 ) = f_;
    
    xrang = wid*[1/3,2/3];
    c_ = code( x >= xrang(1) & x <= xrang(2) & y >= yrang(1) & y <= yrang(2), : );
    f_ = max(c_);
    feas( 1:length(f_), 7 ) = f_;

    xrang = wid*[2/3,1];
    c_ = code( x >= xrang(1) & x <= xrang(2) & y >= yrang(1) & y <= yrang(2), : );
    f_ = max(c_);
    feas( 1:length(f_), 8 ) = f_;










