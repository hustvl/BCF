% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function [Dotau, A, E, f, tfm_matrix, error_sign, UData, VData, XData, YData, A_scale]=tilt_kernel(input_image, mode, center, focus_size, initial_tfm_matrix, para)
% tilt_kernel aligns a sub-image of input_image specified by base_points
% and focus_size to its frontal, with low-column rank subject to some
% linear constraints.
% -------------------------input------------------------------------------
% input_image:  height-by-width real matrix or height*width*3 real matrix but
%               we will only preserve the first channel for the second
%               case.
% mode:         one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 'homography',
%               'homography_notranslation'.
% center:       2-by-1 real vector in row-column co-ordinates, origin of
%               the image.
% focus_size:   1-by-2 real vector in row-column co-ordinates, size of the
%               focus.
% initial_tfm_matrix:   3-by-3 matrix, optional, the initial transform
%                       matrix. If not specified set it to be eye(3);
% para:         parameters for both inner-loop and outer-loop, must be
%               specified!
%
% -------------------------output-----------------------------------------
% Dotau:        real matrix with same size as focus_size, aligned images.
% A:            low-rank part of Dotau;
% E:            sparse-error part of Dotau;
% f:            value of objective-function;
% tfm_matrix:   resulted transform matrix.
% error_sign:   0 or 1, 1 for trival solutions.
%
% -------------------------note-------------------------------------------
% 1. Unless explicitly stated, point follows x-y coordinates and is a
% column vector, and size follows row-column coordinates and is a row
% vector.

disp(['mode=', num2str(mode), ', scale=[', num2str(size(input_image, 1)), '><', num2str(size(input_image, 2))]);
%% parse parameters.
original_image=input_image;
if size(input_image, 3)>1
    input_image=input_image(:, :, 1)*0.299+input_image(:, :, 2)*0.587+input_image(:, :, 3)*0.144;
%     input_image=input_image(:, :, 1);
end
input_image=double(input_image);
% make base_points and focus_size integer, the effect of this operations
% remains to be tested.
center=floor(center);
focus_size=floor(focus_size);
if isempty(para.outer_tol)
    outer_tol=5e-5;
else
    outer_tol=para.outer_tol;
end

if isempty(para.outer_max_iter)
    outer_max_iter=50;
else
    outer_max_iter=para.outer_max_iter;
end

if isempty(para.outer_display_period)
    outer_display_period=1;
else
    outer_display_period=para.outer_display_period;
end

inner_para=[];
if isempty(para.inner_tol)
    inner_para.tol=1e-7;
else
    inner_para.tol=para.inner_tol;
end

if isempty(para.inner_c)
    inner_para.c=1;
else
    inner_para.c=para.inner_c;
end

if isempty(para.inner_mu)
    inner_para.mu=[];
else
    inner_para.mu=para.inner_mu;
end

if isempty(para.inner_display_period)
    inner_para.display_period=100;
else
    inner_para.display_period=para.inner_display_period;
end

if isempty(para.inner_max_iter)
    inner_para.max_iter=inf;
else
    inner_para.max_iter=para.inner_max_iter;
end

if isempty(para.display_result)
    display_result=1;
else
    display_result=para.display_result;
end
%% decide origin of the two axes.
image_size=size(input_image);
image_center=floor(center);
focus_center=zeros(2, 1);
focus_center(1)=floor((1+focus_size(2))/2);
focus_center(2)=floor((1+focus_size(1))/2);
UData=[1-image_center(1) image_size(2)-image_center(1)];
VData=[1-image_center(2) image_size(1)-image_center(2)];
XData=[1-focus_center(1) focus_size(2)-focus_center(1)];
YData=[1-focus_center(2) focus_size(1)-focus_center(2)];
A_scale=1;

Dotau_series={};

%% prepare initial data
input_du=imfilter(input_image, -fspecial('sobel')'/8);
input_dv=imfilter(input_image, -fspecial('sobel')/8);
tfm_matrix=initial_tfm_matrix;
tfm=fliptform(maketform('projective', tfm_matrix'));
Dotau=imtransform(input_image, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
Dotau_series{1}=Dotau;
initial_image=Dotau;
du=imtransform(input_du, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
dv=imtransform(input_dv, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
du= du/norm(Dotau, 'fro')-(sum(sum(Dotau.*du)))/(norm(Dotau, 'fro')^3)*Dotau;
dv= dv/norm(Dotau, 'fro')-(sum(sum(Dotau.*dv)))/(norm(Dotau, 'fro')^3)*Dotau;
A_scale=norm(Dotau, 'fro');
Dotau=Dotau/norm(Dotau, 'fro');

tau=tfm2para(tfm_matrix, XData, YData, mode);
J=jacobi(du, dv, XData, YData, tau, mode);
S=constraints(tau, XData, YData, mode);

outer_round=0;
pre_f=0;
%% begin main loop.
while 1
   outer_round=outer_round+1;
   [A, E, delta_tau, f, error_sign]=inner_IALM_constraints(Dotau, J, S, inner_para);
   if error_sign==1
       return;
   end
   %% display information
   if mod(outer_round, outer_display_period)==0
      disp(['outer_round ',num2str(outer_round), ', f=',num2str(f), ', rank(A)=', num2str(rank(A)), ', ||E||_1=', num2str(sum(sum(abs(E))))]);
   end
   %% update Dotau.
   tau=tau+delta_tau;
   tfm_matrix=para2tfm(tau, XData, YData, mode);
   tfm=fliptform(maketform('projective', tfm_matrix'));
   Dotau=imtransform(input_image, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
   Dotau_series{outer_round+1}=Dotau;
   %% judge convergence
   if outer_round>=outer_max_iter || abs(f-pre_f)<outer_tol
       break;
   end
   %% record data and prepare for the next round.
   pre_f=f;
   du=imtransform(input_du, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
   dv=imtransform(input_dv, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'size', focus_size);
   du= du/norm(Dotau, 'fro')-(sum(sum(Dotau.*du)))/(norm(Dotau, 'fro')^3)*Dotau;
   dv= dv/norm(Dotau, 'fro')-(sum(sum(Dotau.*dv)))/(norm(Dotau, 'fro')^3)*Dotau;
   A_scale=norm(Dotau, 'fro');
   Dotau=Dotau/norm(Dotau, 'fro');
   J=jacobi(du, dv, XData, YData, tau, mode);
   S=constraints(tau, XData, YData, mode);
end

%% display results
if para.display_result==1
   figure(para.figure_no);
   clf;
   subplot(2, 2, 1);
   if max(max(max(initial_image)))~=0
       imshow(initial_image, [], 'DisplayRange', [0 max(max(max(initial_image)))]);
   end
   subplot(2, 2, 2);
   if max(max(max(Dotau)))~=0
       imshow(Dotau, [], 'DisplayRange', [0 max(max(max(Dotau)))]);
   end
   subplot(2, 2, 3);
   if max(max(max(A)))~=0
       imshow(A, [], 'DisplayRange', [0 max(max(A))]);
   end
   subplot(2, 2, 4);
   if max(max(max(E)))~=0
       imshow(E, [], 'DisplayRange', [0 max(max(E))]);
   end
end

if ~isempty(para.save_path)
    if ~exist(para.save_path) || ~isdir(para.save_path)
        mkdir(para.save_path);
    end
    imwrite(uint8(original_image), fullfile(para.save_path, 'input.jpg'), 'JPEG');
    imwrite(uint8(Dotau), fullfile(para.save_path, 'Dotau.jpg'), 'JPEG');
    imwrite(uint8(A*A_scale), fullfile(para.save_path, 'A.jpg'), 'JPEG');
    imwrite(uint8(abs(E)/max(max(abs(E)))*255), fullfile(para.save_path, 'E.jpg'), 'JPEG');
    for i=1:length(Dotau_series)
        imwrite(uint8(abs(Dotau_series{i})/max(max(abs(Dotau_series{i})))*255), fullfile(para.save_path, ['Dotau_', num2str(i), '.jpg']), 'JPEG');
    end
end