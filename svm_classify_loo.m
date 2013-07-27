function svm_classify_loo(fea_dir, c_svm)
% -------------------------------------------------------------------------
% evaluate the performance of the image feature using linear SVM
% we used Liblinear package in this example code

fdatabase = retr_database_dir(fea_dir, '*.mat');

fprintf('\nLeave one out testing...\n');

len = length( fdatabase.path );
feas = zeros( len, 0, 'single' );
labels = zeros( len, 1, 'single' );

for i = 1:len
    fprintf('Loading pre computed features %d of %d.\n', i, len);
    
    fpath = fdatabase.path{i};
    load(fpath, 'fea', 'label');
    
    feas(i, 1:length(fea) ) = full(fea);
    labels(i) = label;
end

options = ['-s 4 -c ' num2str(c_svm)];
accuracy = zeros( len, 1 );

n_class = length(unique( labels ));
i_all_1 = cell(n_class, 1);
i_all   = cell(n_class, 1);
for i = 1:n_class
    i_c = find( labels == i );
    i_all{i} = i_c;
    rid = randperm( length(i_c) );
    i_all_1{i} = i_c( rid(1:end-1) );
end

for i = 1:len
    fprintf( 'Testing %d of %d.\n', i, len );
    
    feate = sparse(double(feas(i, : )));
    labte = sparse(double(labels(i)));
    
    indtr = i_all_1;
    indtr{labte} = setdiff( i_all{labte}, i );
    indtr = cell2mat( indtr );
    
    model = train( labels(indtr), feas(indtr,:) , options);
    [~, accuracy(i)] = predict( labte, feate, model );
end

Ravg = mean(accuracy/100);                  % average recognition rate
Rstd = std(accuracy/100);                   % standard deviation of the recognition rate

fprintf('\n===============================================\n');
fprintf('Average classification accuracy: %f\n', Ravg);
fprintf('Standard deviation: %f\n', Rstd);    
fprintf('===============================================\n');
