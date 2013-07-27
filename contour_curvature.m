
function cs = contour_curvature( cont )

    ispl = 1 : 10 : size(cont, 1);
    cs = zeros( length(ispl), 1 );
    
    for i = 1:length( ispl )
        cs(i) = curvature( cont, ispl(i), 5 );
        if isnan( cs(i) )
            disp('sss');
        end
    end