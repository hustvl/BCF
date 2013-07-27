function [segx, segy]=generate_segments_2(c, maxvalue, N, nn)

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
n_kp = size( kp, 1 );

n2 = dist2(kp, c);
[~, i_kp] = min(n2'); 

n_cf = (n_kp-1)*n_kp/2;
segx = zeros( n_cf, nn );
segy = zeros( n_cf, nn );

s = 1;
for i = 1:n_kp
    for j = i+1:n_kp
        cf = c( i_kp(i) : i_kp(j), : );
        % [XIs,YIs]	= uniform_interp(cf(:,1),cf(:,2),nn-1);
        
        segx(s, : ) = [ cf(1,1); XIs ];
        segy(s, : ) = [ cf(1,2); YIs ];
        s = s + 1;
    end
end



