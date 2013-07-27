% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

% hello_affine.m will show you how to run TILT affine on an image, step by
% step.

clc;
close all;
task_name='TILT_AFFINE';

%% load the image
[img_name, img_path]=uigetfile('*.*');
img=imread(fullfile(img_path, img_name));

%% get the top-left and bottom-right corner of the rectangle window where
%% we perform TILT.
% by default, the first point should be top-left one and the second should
% be bottom-right. Don't mess up the order......
figure(1);
imshow(img);
hold on;
initial_points=zeros(2, 2);
for i=1:2
    initial_points(:, i)=ginput(1)';
    plot(initial_points(1, i), initial_points(2, i), 'rx');
    hold on;
end

%% Run TILT:
% Suppose the input image is an m0*n0 matrix, and the transformed window is
% m1*n1. First to establish the mapping we have to specify the coordinates
% that the two image domains are in. Naively for a h*w image, we should
% place the top-left corner at (0, 0) and bottom right corner at (w-1,
% m-1), but here for convenience we choose a different coordinate, which is
% specified below. Suppose the input image lies in the uv coordinate(u for
% width and v for height), and the transformed window lies in xy
% coordinate. Then we put the two corners of the input image at point
% (UData(1), VData(1))^T and (UData(2), VData(2))^T, and similarly we put
% the transformed window at (XData(1), YData(1))^T and (XData(2),
% YData(2))^T. 
%
% After specifying the coordinate, we have to describe the mapping. For
% some reason, we use tfm_matrix to describe not the forward map from the
% input image to the transformed window, but the inverse one, i.e., from
% the transformed window to the input image. So when you want to use
% imtransform, please remember to flipt the tfm you get......
tic

% Here the first three parameters are required in all cases. To tune
% parameters, you can set value of the parameter that you want to tune in
% the form ...,'param_name', para_value,...
% Normally the default parameters acts pretty well.

% Introduce in some important parameters:
% 1.BLUR_SIGMA, and BLUR_NEIGHBOR: Indeed these two parameters controls how
% large you blur an image. If the blur is not enough, then the convergence
% range of TILT may be quite small. On the contrary if you blur too much,
% then all the textures may be blurred out then you have no hope to recover
% the structure. If something is wrong with TILT, first please tune this
% parameter.
%
% 2.PYRAMID_MAXLEVEL, normally if the texture is rich enough, running the
% TILT in the lowest-resolution level should already give satisfactory
% result. TILT on the higher-level is much slower and may not improves too
% much.
%
% 3.SAVE_PATH, to make it easy to use, TILT can save some intermediate
% results onto the disk to give users an idea about what's going on. Please
% speicyfing the full path(absolute path), not the relative path of the
% direction that you want to save. If you do not specify this, by default,
% TILT will not save these intermediate results.
%
% 4.DISPLAY_INTER, 1 if you want to show the intermediate results, 0 if
% not. If you just run a large number of TILT, perhaps setting this to
% 0(disabling intermediate displaying) would accelerate......
[Dotau, A, E, f, tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
    TILT(img, 'AFFINE', initial_points, 'SAVE_PATH', fullfile(cd, task_name),...
    'BLUR_SIGMA', 2, 'BLUR_NEIGHBOR', 2,...
    'PYRAMID_MAXLEVEL', 1, 'DISPLAY_INTER', 1);
toc
%% examples on how to use the result of TILT:
% 1. transform the input image to the focus window.
tfm=fliptform(maketform('projective', tfm_matrix'));
Dotau=imtransform(img, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'Size', focus_size);
figure(3);
imshow(Dotau);
% 2. rectify the input image.XData_whole, YData_whole
% specify the location of the whole rectified image in x-y coordinate.
[rectify_image, XData_whole, YData_whole]=imtransform(img, tfm, 'bilinear', 'UData', UData, 'VData', VData);
figure(4);
imshow(rectify_image);
% 3. adjust tfm_matrix to standard coordinate, 
standard_tfm_matrix=[1 0 1-UData(1); 0 1 1-VData(1); 0 0 1]*tfm_matrix*[1 0 XData(1)-1; 0 1 YData(1)-1; 0 0 1]; % indeed just compose the tfm_matrix with coordinate swtich.
standard_tfm=fliptform(maketform('projective', standard_tfm_matrix'));
Dotau=imtransform(img, standard_tfm, 'bilinear', 'UData', [1 size(img, 2)], 'VData', [1 size(img, 1)],  'XData', [1 focus_size(2)], 'YData', [1 focus_size(1)], 'Size', focus_size);
figure(5);
imshow(Dotau);