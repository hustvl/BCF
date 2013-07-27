function feat = shape_signature(pnts)


    x0 = pnts(:, 1);
    y0 = pnts(:, 2);
    
    x = x0 - x0(1);
    y = y0 - y0(1);
    
    s = 1 / sqrt(x(end)^2 + y(end)^2);
    x = x * s;
    y = y * s;
        
    x2 = x(end) * x + y(end) * y;
    y2 = -y(end) * x + x(end) * y;
    
    feat = [x2; y2]';
        
   