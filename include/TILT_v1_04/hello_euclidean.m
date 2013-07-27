% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

% hello_euclidean.m will show you how to run TILT with euclidean mode on an image, step by
% step.

clc;
close all;
task_name='TILT_EUCLIDEAN';

%% load the image
clc;
clear;
close all;
b_reload=1;
if b_reload
    [img_name, img_path]=uigetfile('*.*');
    img=imread(fullfile(img_path, img_name));
    
    %%% synthetic test on TILT performed on large images.

    %%% add synthetic rotation to the image.
    theta=pi/32;
    if theta~=0
        UData=[-size(img, 2)/2 size(img, 2)/2];
        VData=[-size(img, 1)/2 size(img, 1)/2];
        XData=UData;
        YData=VData;
        tfm_matrix=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        img=imtransform(img, maketform('projective', tfm_matrix'), 'bilinear', ...
            'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'Size', [size(img, 1) size(img, 2)]);
    end

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
%%% code for validating the effectiveness of the optimization,
%%% branch-and-band here is closed.
% tic
% [Dotau, A, E, f, tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
%     TILT(img, 'EUCLIDEAN', initial_points, 'BRANCH', 0, 'SAVE_PATH', [], 'DISPLAY_INTER', 1);
% toc

%%% code for customers.
tic
[Dotau, A, E, f, tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
    TILT(img, 'EUCLIDEAN', initial_points, 'SAVE_PATH', [], 'DISPLAY_INTER', 1);
toc