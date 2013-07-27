clear 

addpath('include/llc/');
addpath('include');
addpath('include/liblinear-1.93/matlab'); 
addpath('include/idsc_distribute/common_innerdist/');
addpath('include/vl_toolbox');
vl_setup;
addpath(genpath('include/edge_re_linking/'))

img_dir = '../101_ObjectCategories/';
fea_dir = 'data/101_feas';
codebook_path = 'data/101_codebook.mat';
fmt = '*.bmp';


para.n_shapesamp = 2000;
para.n_contsamp = 50;
para.max_curvature = 0.5;
para.n_pntsamp = 50;
para.k_sc = 1024;
para.knn = 5;
para.t_edge = 0.1;
para.min_len_cf = 15;	% filter out the very short edges


C = 10;
tr_num = 30;                % number of training images
nRounds = 10;
pyramid = [1,2,3,4];          %[1,2,3,4];

extr_cf_101(img_dir, '*.bmp', para);
learn_codebook(img_dir, codebook_path, para);

encode_cf(img_dir, codebook_path, para);
encode_cf_fake(img_dir, codebook_path, para);

pyramid_mmpooling(img_dir, fea_dir, pyramid);
svm_classify(fea_dir, tr_num, C, nRounds, inf);



