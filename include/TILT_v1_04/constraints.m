% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function S=constraints(tau, XData, YData, mode)
% constraints() will get the linearize constraints of tau according to
% mode.
% -----------------------------input--------------------------------------
% tau:          p-by-1 real vector.
% mode:         one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 'homography',
%               'homography_notranslation'.
% ----------------------------output--------------------------------------
% S:            linearized constraints on tau.
switch lower(mode)
    case 'euclidean'
        S=zeros(2, 1);
    case 'affine'
        S=zeros(2, 4);
        vec1=[tau(1);tau(3)];
        vec2=[tau(2);tau(4)];
        V=sqrt(norm(vec1)^2*norm(vec2)^2-(vec1'*vec2)^2);
        S(1, 1)=(tau(1)*norm(vec2)^2-vec1'*vec2*tau(2))/V;
        S(1, 2)=(tau(2)*norm(vec1)^2-vec1'*vec2*tau(1))/V;
        S(1, 3)=(tau(3)*norm(vec2)^2-vec1'*vec2*tau(4))/V;
        S(1, 4)=(tau(4)*norm(vec1)^2-vec1'*vec2*tau(3))/V;
        S(2, 1)=2*tau(1);
        S(2, 2)=-2*tau(2);
        S(2, 3)=2*tau(3);
        S(2, 4)=-2*tau(4);
    case 'homography'
        S=zeros(1, 8);
        temp=reshape(tau, 2, 4);
        X=temp(1, :);
        Y=temp(2, :);
        e1=[X(3)-X(1); Y(3)-Y(1)];
        e2=[X(4)-X(2); Y(4)-Y(2)];
        norm_e1=e1'*e1;
        norm_e2=e2'*e2;
        e1e2=e1'*e2;
        N=2*sqrt(norm_e1*norm_e2-e1e2^2);
        S(1, 1)=1/N*(2*(X(3)-X(1))*(-1)*norm_e2-2*e1e2*(-1)*(X(4)-X(2)));
        S(1, 2)=1/N*(2*(Y(3)-Y(1))*(-1)*norm_e2-2*e1e2*(-1)*(Y(4)-Y(2)));
        S(1, 3)=1/N*(2*(X(4)-X(2))*(-1)*norm_e1-2*e1e2*(-1)*(X(3)-X(1)));
        S(1, 4)=1/N*(2*(Y(4)-Y(2))*(-1)*norm_e1-2*e1e2*(-1)*(Y(3)-Y(1)));
        S(1, 5)=1/N*(2*(X(3)-X(1))*norm_e2-2*e1e2*(X(4)-X(2)));
        S(1, 6)=1/N*(2*(Y(3)-Y(1))*norm_e2-2*e1e2*(Y(4)-Y(2)));
        S(1, 7)=1/N*(2*(X(4)-X(2))*norm_e1-2*e1e2*(X(3)-X(1)));
        S(1, 8)=1/N*(2*(Y(4)-Y(2))*norm_e1-2*e1e2*(Y(3)-Y(1)));
end