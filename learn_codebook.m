function learn_codebook(img_dir, codebook_path, para)

    database = retr_database_dir(img_dir, '*_feat.mat');
    N = length(database.path);

    ri = randperm(N);
    M = N;
    feats_sc = cell( 1, M );
    for n = 1:M
        fprintf('Loading shape context features: %d of %d\n', n, M );
        load(database.path{ri(n)}, 'feat_sc' );    

        nmax = 800;
        if size(feat_sc, 2) > nmax
            ri2 = randperm( size(feat_sc,2) );
            feat_sc = feat_sc(:, ri2(1:nmax) );
        end

        feats_sc{n} = feat_sc;
    end
    feats_sc = single(cell2mat(feats_sc));

    [dict_sc, assign] = vl_kmeans(feats_sc, para.k_sc, 'verbose', 'distance', 'l2', 'algorithm', 'elkan');
    dict.sc = dict_sc;

    save(codebook_path, 'dict' );







