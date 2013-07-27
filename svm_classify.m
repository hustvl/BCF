function [Ravg, Rstd, acct] = svm_classify(fea_dir, tr_num, c_svm, nRounds, mem_block)
% -------------------------------------------------------------------------
% evaluate the performance of the image feature using linear SVM
% we used Liblinear package in this example code

fdatabase = retr_database_dir(fea_dir, '*.mat');


fprintf('\n Testing...\n');
clabel = unique(fdatabase.label);
nclass = length(clabel);
accuracy = zeros(nRounds, 1);

acct = [];
wpath = {};

as = zeros( nclass, nRounds );

for ii = 1:nRounds,
    fprintf('Round: %d...\n', ii);
    tr_idx = [];
    ts_idx = [];
    
    for jj = 1:nclass,
        idx_label = find(fdatabase.label == clabel(jj));
        num = length(idx_label);
        
        idx_rand = randperm(num);
      
        tr_idx = [tr_idx; idx_label(idx_rand(1:tr_num))];
        ts_idx = [ts_idx; idx_label(idx_rand(tr_num+1:end))];
    end
    
    fprintf('Training number: %d\n', length(tr_idx));
    
    % load the training features 
    tr_fea = zeros(length(tr_idx), 0, 'single');
    tr_label = zeros(length(tr_idx), 1, 'single');
    
    fprintf('Loading training data....\n');
    for jj = 1:length(tr_idx),
        fpath = fdatabase.path{tr_idx(jj)};
        load(fpath, 'fea', 'label');
       
        tr_fea(jj, 1:length(fea) ) = full(fea);
        tr_label(jj) = label;
    end
    
    fprintf('svm training....\n');
    tic
    options = ['-s 4 -c ' num2str(c_svm)];
    model = train(tr_label, tr_fea, options);
    clear tr_fea;
    toc

    fprintf('Testing number:%d\n', length(ts_idx));
    % load the testing features
    ts_num = length(ts_idx);
    ts_label = [];
    
    if ts_num < mem_block,
        % load the testing features directly into memory for testing
        ts_fea = zeros(length(ts_idx), 0);
        ts_label = zeros(length(ts_idx), 1);

        for jj = 1:length(ts_idx),
            fpath = fdatabase.path{ts_idx(jj)};
            load(fpath, 'fea', 'label');
            ts_fea(jj, 1:length(fea)) = full(fea);
            ts_label(jj) = label;
        end

        [C] = predict(ts_label, sparse(ts_fea), model);
    else
        % load the testing features block by block
        num_block = floor(ts_num/mem_block);
        rem_fea = rem(ts_num, mem_block);
        
        curr_ts_fea = zeros(mem_block, dFea);
        curr_ts_label = zeros(mem_block, 1);
        
        C = [];
        
        for jj = 1:num_block,
            block_idx = (jj-1)*mem_block + (1:mem_block);
            curr_idx = ts_idx(block_idx); 
            
            % load the current block of features
            for kk = 1:mem_block,
                fpath = fdatabase.path{curr_idx(kk)};
                load(fpath, 'fea', 'label');
                curr_ts_fea(kk, :) = full(fea);
                curr_ts_label(kk) = label;
            end    
            
            % test the current block features
            ts_label = [ts_label; curr_ts_label];
            [curr_C] = predict(curr_ts_label, sparse(curr_ts_fea), model);
            C = [C; curr_C];
        end
        
        curr_ts_fea = zeros(rem_fea, dFea);
        curr_ts_label = zeros(rem_fea, 1);
        curr_idx = ts_idx(num_block*mem_block + (1:rem_fea));
        
        for kk = 1:rem_fea,
           fpath = fdatabase.path{curr_idx(kk)};
           load(fpath, 'fea', 'label');
           curr_ts_fea(kk, :) = fea;
           curr_ts_label(kk) = label;
        end  
        
        ts_label = [ts_label; curr_ts_label];
        [curr_C] = predict(curr_ts_label, sparse(curr_ts_fea), model); 
        C = [C; curr_C];        
    end
    
    acc = C == ts_label;

    accuracy(ii) = mean(acc); 
    fprintf('Classification accuracy for round %d: %f\n', ii, accuracy(ii));
    
    acct = [acct, acc];
    
    % find the mis-classified shapes
    % wrongs = find( C ~= ts_label );
    % for w = 1:length(wrongs)
    %     wpath{ length(wpath)+1 } = fdatabase.path{ wrongs(w) };
    % end
    
    for ic = 1:nclass
        as( ic, ii ) = mean( acc(ts_label == ic) );
    end
    
end

Ravg = mean(accuracy);                  % average recognition rate
Rstd = std(accuracy);                   % standard deviation of the recognition rate

fprintf('\n===============================================\n');
fprintf('Average classification accuracy: %f\n', Ravg);
fprintf('Standard deviation: %f\n', Rstd);    
fprintf('===============================================\n\n');


% whist = zeros(1, nclass);
% for i = 1:length(wpath)
%     [~, fn] = fileparts( wpath{i} );
%     whist( str2num(fn(1:2)) ) = whist( str2num(fn(1:2)) ) + 1;
%     % fprintf('%s\n', wpath{i});
% end
% whist = whist / sum(whist);
% figure, bar(whist);
