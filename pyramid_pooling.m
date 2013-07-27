function pyramid_pooling(img_dir, fea_dir, pyramid)

database = retr_database_dir(img_dir, '*_code.mat');


for  n = 1:length( database.path )
    
    label = database.label(n);
    path = database.path{n};
    load(path);
    
    fprintf('pyramid feature context %d of %d\n', n , length(database.path));
    feas = zeros( 0, sum(pyramid.^2) );
    counter = 0;
    
    for p = 1:length(pyramid)
        
        for i = 1:pyramid(p)
            for j = 1:pyramid(p)
                
                yrang = fea_rbf.hgt*[i-1,i]/pyramid(p);
                xrang = fea_rbf.wid*[j-1,j]/pyramid(p);
                
                
                c_ = fea_rbf.code( fea_rbf.x >= xrang(1) & fea_rbf.x <= xrang(2) & fea_rbf.y >= yrang(1) & fea_rbf.y <= yrang(2), : );
                f_ = max(c_);
                
                counter = counter + 1;
                feas( 1:length(f_), counter ) = f_;
            end
        end
    
    end

    fea = feas(:);
    fea = fea / sqrt(sum(fea.^2));
    
    [~, fname] = fileparts(path);
    feapath = fullfile(fea_dir, num2str(label), [fname(1:end-5) '.mat']);
    if ~isdir(fullfile(fea_dir, num2str(label))),
        mkdir(fullfile(fea_dir, num2str(label)));
    end  
    
    save(feapath, 'fea', 'label');

    
end