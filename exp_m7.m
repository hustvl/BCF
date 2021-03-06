
% Experiment on Mpeg-7 dataset

clear; clc;

addpath('include/llc/');
addpath('include');
addpath('include/liblinear-1.7-single/matlab'); 
addpath('include/idsc_distribute/common_innerdist/');
addpath('include/vl_toolbox');
vl_setup;

img_dir = 'data/mpeg71';
fea_dir = 'data/mpeg71_feas';
codebook_path = 'data/mpeg71_codebook.mat';
fmt = '*.bmp';


para.n_shapesamp = 2000;
para.n_contsamp = 50;
para.max_curvature = 0.5;
para.n_pntsamp = 100;
para.k_sc = 1500;
para.knn = 5;

C = 10;
tr_num = 10;                % number of training images
nRounds = 10;
pyramid = [1, 2, 4];          %[1,2,3,4];

normalize_shape(img_dir, '*.bmp');
extr_cf(img_dir, '*.png', para);
learn_codebook(img_dir, codebook_path, para);
encode_cf(img_dir, codebook_path, para);
pyramid_mmpooling(img_dir, fea_dir, pyramid);

svm_classify(fea_dir, tr_num, C, nRounds, inf);		% half testing
svm_classify_loo(fea_dir, C);						% leave one out 

