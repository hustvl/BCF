% globalcurv
% Longin Jan Latecki January 2007
%
% This function computes the curvature for each point in 
% a X-by-2 matrix of points
%
% Parameters:
%    1. P - the X-by-2 matrix of points
%    2. nsize - the number of neighbor points to use on each side of a
%       point
%    3. SEGMENT - a boolean variable that is set to 1 if the points
%       represent a segment, and set to 0 if the points represent a closed
%       curve
%
% Return values:
%    1. gc - the vector of the curvature values
%       Note: if the points represent a segment, then the size of gc will
%       be 2*nsize smaller than number of points in P. If the points
%       represent a closed curve, then gc will be the same size as the
%       number of points in P
%
% Note: if the points represent a segment, then the curvature is not
% calculated for the first nsize points nor the last nsize points

% Example:
%    input example: P is 100x2 vector of points coordinates on the plane
%    nsize=5, giving for each point the neighborhood of 11 points

function gc=globalcurv(P, nsize, SEGMENT);
    len = size(P,1);
    gc = zeros(1, len);
    if ( SEGMENT == 1 ) %The points represent a segment
        for i = nsize+1:len-nsize
            gc(i-nsize) = curvature( P, i, nsize );
        end
        gc(1:nsize) = gc(1);
        gc(len - nsize + 1: len) = gc(len - nsize);
    else
        for i = 1:len
            gc(i) = curvature( P, i, nsize );
        end
    end
return