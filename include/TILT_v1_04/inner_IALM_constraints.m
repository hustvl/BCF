% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function [A, E, delta_tau, f, error_sign]=inner_IALM_constraints(Dotau, J, S_J, inner_para)
% inner_IALM_constraints will solve the programming:
% min ||A||_*+lambda*||E||_1   s.t.     Dotau+J*delta_tau=A+E and
%                                       S*delta_tau=0
% via the Inexact ALM method.
% 
% ---------------------------------input----------------------------------
% Dotau:            m-by-n image matrix.
% J:                m-by-n-by-p tensor.
% S_J:                c-by-p matrix.
% inner_para:       parameters.
% --------------------------------output----------------------------------
% A:                m-by-n matrix, low-rank part.
% E:                m-by-n matrix, error part.
% delta_tau:        step of tau.
% f:                objective funtion value.

%% Parse Inputs.
try
    if isempty(inner_para.tol)
        tol=1e-7;
    else
        tol=inner_para.tol;
    end

    if isempty(inner_para.c)
        c=1;
    else
        c=inner_para.c;
    end

    if isempty(inner_para.mu)
        mu=1.25/norm(Dotau);
    else
        mu=inner_para.mu;
    end

    if isempty(inner_para.display_period)
        display_period=100;
    else
        display_period=inner_para.display_period;
    end

    if isempty(inner_para.max_iter)
        max_iter=inf;
    else
        max_iter=inner_para.max_iter;
    end

    %% prepare data
    [m n]=size(Dotau);
    E=zeros(m, n);
    A=zeros(m, n);
    p=size(J, 3);
    delta_tau=zeros(p, 1);

    J_vec=reshape(J, m*n, p);
    Jo=J_vec;
    J_vec=[J_vec; S_J];
    pinv_J_vec=pinv(J_vec);
    inner_round=0;
    rho=1.25;
    lambda=c/sqrt(m);

    Y_1=Dotau;
    norm_two=norm(Y_1, 2);
    norm_inf=norm(Y_1(:), inf)/lambda;
    dual_norm=max(norm_two, norm_inf);
    Y_1=Y_1/dual_norm;
    Y_2=zeros(size(S_J, 1), 1);
    d_norm=norm(Dotau, 'fro');
    error_sign=0;
    first_f=sum(svd(Dotau));
catch
    error_sign=1;
    A=Dotau;
    E=zeros(size(Dotau));
    delta_tau=zeros(size(J, 3), 1);
    f=inf;
    return;
end


%% begin main loop
while 1
    try
        inner_round=inner_round+1;
        temp_0=Dotau+reshape(Jo*delta_tau, m, n)+Y_1/mu;
        temp_1=temp_0-E;
        [U S V]=svd(temp_1, 'econ');
        A=U*((S>1/mu).*(S-1/mu))*V';
        temp_2=temp_0-A;
        %E_new=(temp_2>lambda/mu).*(temp_2-lambda/mu)+(temp_2<-lambda/mu).*(temp_2+lambda/mu);
		E=(temp_2>lambda/mu).*(temp_2-lambda/mu)+(temp_2<-lambda/mu).*(temp_2+lambda/mu);
        f=sum(sum(abs((S>1/mu).*(S-1/mu))))+lambda*sum(sum(abs(E)));
        temp_3=A+E-Dotau-Y_1/mu;
        temp_3=reshape(temp_3, m*n, 1);
        temp_3=[temp_3; -Y_2/mu];
        delta_tau=pinv_J_vec*temp_3;
        derivative_Y_1=Dotau-A-E+reshape(Jo*delta_tau, m, n);
        derivative_Y_2=S_J*delta_tau;
        Y_1=Y_1+derivative_Y_1*mu;
        Y_2=Y_2+derivative_Y_2*mu;
    catch
        error_sign=1;
        A=Dotau;
        E=zeros(size(Dotau));
        delta_tau=zeros(size(J, 3), 1);
        f=first_f;
        return;
    end
    %% judge error
    if f<first_f/3
        error_sign=1;
        A=Dotau;
        E=0;
        f=first_f;
        delta_tau=zeros(size(J, 3), 1);
        return;
    end
    %% display
    stop_criterion=sqrt(norm(derivative_Y_1, 'fro')^2+norm(derivative_Y_2, 2)^2)/d_norm;
    if stop_criterion<tol || inner_round >max_iter
        break;
    end
    if mod(inner_round, display_period)==0
        disp(['        inner round ', num2str(inner_round), ' stop_criterion=',num2str(stop_criterion), 'rank(A)=',num2str(rank(A_new)), '||E||_1=',num2str(sum(sum(abs(E_new))))]);
    end

    %% update
    mu=mu*rho;
end
