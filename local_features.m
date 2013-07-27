function [feat, xy, cfs] = local_features(c, maxvalue, N, nn)

%-------------------------------------------------------
%[SegmentX, SegmentY,NO]=GenSegmentsNew(a,b,maxvalue,nn)
%This function is used to generate all the segments 
%vectors of the input contour
%a and b are the input contour sequence
% maxvalue is the stop condition of DCE, usually 1~1.5
% nn is the sample points' number on each segment, in super's method,n=25

%SegmentX,SegmentY denotes all the coordinates of all the segments of input contour
%NO denotes the number of segments of input contour 
%-------------------------------------------------------


kp = evolution(c, N, maxvalue, 0, 0, 1);

n2 = dist2(kp, c);
[~, i_kp] = min(n2'); 

% remove some very short segments
% d = i_kp - [-inf, i_kp(1:end-1)];
% i_kp(d<50) = [];


n_kp = length(i_kp);

n_cf = (n_kp-1)*n_kp/2;

feat = zeros( 0, n_cf );
xy = zeros( n_cf, 2 );
cfs = cell( n_cf, 1 );

s = 1;
for i = 1:n_kp
    for j = i+1:n_kp
        cf = c( i_kp(i) : i_kp(j), : );
        
        len_ = size( cf, 1 );
        ii_ = round( 1 : (len_-1)/(nn-1) : len_ );
        cf = cf(ii_,:);
        
        sc = shape_context( cf );
        sc = sc(:);
        sc = sc / (sum(sc)+eps);
        feat( 1:length(sc), s ) = sc;
        
        xy( s, : ) = cf(round( nn/2 ), : );
        
        cfs{s} = cf;
        
        s = s + 1;
    end
end



