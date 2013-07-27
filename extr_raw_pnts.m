function pnts = extr_raw_pnts(c, maxvalue, N, nn)

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

n_kp = length(i_kp);

n_cf = (n_kp-1)*n_kp + 1;

pnts = cell( n_cf, 1 );

s = 1;
for i = 1:n_kp
    for j = 1:n_kp
        
        if i == j 
            continue;
        end
        
        if i < j
            cf = c( i_kp(i) : i_kp(j), : );
        end
        
        if i > j
            cf = [ c( i_kp(i) : end, : ); c( 1 : i_kp(j), : ) ];
        end
        
        pnts{s} = sample_contour(cf, nn);
        s = s + 1;
    end
end

pnts{s} = sample_contour(c, nn);

function cf = sample_contour( cf, nn )
len_ = size( cf, 1 );
ii_ = round( 1 : (len_-1)/(nn-1) : len_ );
cf = cf(ii_,:);       
        

