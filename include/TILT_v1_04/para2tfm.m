% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function tfm_matrix=para2tfm(tau, XData, YData, mode)
% para2tfm will turn tau to tfm_matrix according to mode.
% ----------------------------input---------------------------------------
% tau:      p-by-1 vector
% mode:     one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 'homography',
%           'homography_notranslation'
% ----------------------------output--------------------------------------
% tfm_matrix:   3-by-3 transform matrix.
tfm_matrix=eye(3);
switch lower(mode)
    case 'euclidean'
        tfm_matrix(1:2, 1:2)=[cos(tau(1)) -sin(tau(1));
                              sin(tau(1)) cos(tau(1))];
    case 'affine'
        tfm_matrix(1, 1:2)=tau(1:2)';
        tfm_matrix(2, 1:2)=tau(3:4)';
    case 'homography'
        X=[XData(1) XData(2) XData(2) XData(1)];
        Y=[YData(1) YData(1) YData(2) YData(2)];
        temp=reshape(tau, 2, 4);
        U=temp(1, :);
        V=temp(2, :);
        A=zeros(8, 8);
        b=zeros(8, 1);
        insert_A=zeros(2, 8);
        insert_b=zeros(2, 1);
        for i=1:4
            insert_A(1, :)=[0 0 0 -X(i) -Y(i) -1 V(i)*X(i) V(i)*Y(i)];
            insert_b(1)=-V(i);
            insert_A(2, :)=[X(i) Y(i) 1 0 0 0 -U(i)*X(i) -U(i)*Y(i)];
            insert_b(2)=U(i);
            A(2*i-1:2*i, :)=insert_A;
            b(2*i-1:2*i)=insert_b;
        end
        solution=A\b;
        tfm_matrix=reshape([solution; 1], 3, 3)';
end
