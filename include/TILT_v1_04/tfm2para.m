% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function tau=tfm2para(tfm_matrix, XData, YData, mode)
% tfm2para will transpose tfm_matrix to its corresponding parameter.
% -------------------------input------------------------------------------
% tfm_matrix:       3-by-3 matrix.
% mode:             one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 
%                   'homography', 'homography_notranslation'
% -------------------------output-----------------------------------------
% tau:              p-by-1 real vector.
switch lower(mode)
    case 'euclidean'
        tau=acos(tfm_matrix(1, 1));
        if tfm_matrix(2, 1)<0
            tau=-tau;
        end
    case 'affine'
        tau=zeros(4, 1);
        tau(1:2)=tfm_matrix(1, 1:2)';
        tau(3:4)=tfm_matrix(2, 1:2)';
    case 'homography'
        X=[XData(1) XData(2) XData(2) XData(1)];
        Y=[YData(1) YData(1) YData(2) YData(2)];
        pt=[X; Y; ones(1, 4)];
        tfm_pt=tfm_matrix*pt;
        tfm_pt(1, :)=tfm_pt(1, :)./tfm_pt(3, :);
        tfm_pt(2, :)=tfm_pt(2, :)./tfm_pt(3, :);
        tau=reshape(tfm_pt(1:2, :), 8, 1);
end