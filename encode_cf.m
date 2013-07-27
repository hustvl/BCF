
function encode_cf(img_dir, codebook_path, para)
    
    database = retr_database_dir(img_dir, '*_feat.mat');
    load(codebook_path);
    dict = dict;
    
    parfor n = 1:length(database.path)
        fprintf('Encoding: %d of %d\n', n, length(database.path));
        func_encode_cf(database.path{n}, dict, para ); 
    end

function func_encode_cf(im_name, dict, para )
    

    load(im_name, 'feat_sc', 'xy', 'sz' );
    
    code_sc = LLC_coding_appr(dict.sc', feat_sc', para.knn)';
    code = code_sc;

    save([im_name(1:end-9), '_code.mat' ], 'code', 'xy', 'sz' );
                                


 
	
	
	
	
	
	