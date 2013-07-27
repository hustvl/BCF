% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function output_pt=transform_point(input_pt, tfm_matrix)
% transform_point will transform input_pt according to tfm_matrix
% ------------------------input-------------------------------------------
% input_pt:         1-by-2 or 2-by-1 real vector
% tfm_matrix:       3-by-3 real matrix.
% ------------------------output------------------------------------------
% output_pt:        a vector of the same size as input_pt.

if size(input_pt, 1)==1
    input_pt=input_pt';
    b_row=1;
else
    b_row=0;
end
input_pt=[input_pt;1];
output_pt=tfm_matrix*input_pt;
output_pt=output_pt/output_pt(3);
output_pt=output_pt(1:2, 1);
if b_row==1
    output_pt=output';
end