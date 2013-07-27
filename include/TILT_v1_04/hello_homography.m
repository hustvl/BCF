% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

% hello_homography.m will show you how to run TILT homography on an image, step by
% step.

clc;
close all;
task_name='TILT_HOMOGRAPHY';

%% load the image
clc;
clear;
close all;
b_reload=1;
if b_reload
    [img_name, img_path]=uigetfile('*.*');
    img=imread(fullfile(img_path, img_name));
    
    %%% synthetic test on TILT performed on large images.
%     img=imresize(img, [size(img, 1)*4, size(img, 2)*4]);

    %% get the top-left and bottom-right corner of the rectangle window where
    %% we perform TILT.
    % by default, the first point should be top-left one and the second should
    % be bottom-right. Don't mess up the order......
    figure(1);
    imshow(img);
    hold on
    initial_points=zeros(2, 2);
    for i=1:2
        initial_points(:, i)=ginput(1)';
        plot(initial_points(1, i), initial_points(2, i), 'rx');
        hold on;
    end
    save('all.mat', 'img', 'initial_points');
end
load all;
%% Run TILT:
tic
[Dotau, A, E, f, tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
    TILT(img, 'HOMOGRAPHY', initial_points, 'SAVE_PATH', [], 'DISPLAY_INTER', 1);
toc