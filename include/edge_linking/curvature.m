% curvature
% Longin Latecki January 2007
%
% This is a function that computes the curvature of a point P(pIndex)
%
% Parameters:
%    1. P - the X-by-2 matrix of points
%    2. index - the index of the point that is currently being worked on
%    3. nbsize - the number of neighbor points to use on each side of the
%       point
%
% Return values:
%    1. c - the curvature calculated for the particular point
%
function c=curvature(P,index,nbsize)
% computes curvature of point P(pIndex)

RADIUS=nbsize;
if size(P,1)<2*RADIUS+1
    error('not enough points');
end

if index<RADIUS+1
    P=[P(end-RADIUS,:); P];
    index=index+RADIUS;
end
if index+RADIUS > size(P,1)
    P=[P;P(1:RADIUS,:)];
end

nb=P(index-RADIUS:index+RADIUS,:);
center=mean(nb);
[U,S,V] = svd(cov(nb));
lam1=S(1,1);

v=center-P(index,:);
dirv=P(index+1,:)-P(index,:);

signv=sign(det([dirv;v]));

c=signv*norm(v)/sqrt(lam1);

