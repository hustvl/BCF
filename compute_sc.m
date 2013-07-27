
function sc_feat = compute_sc( points, nvert, points_num, bins_num, r_inner, r_outer)
% suppose the points are in order

if nargin <= 4
    r_inner=1/8;
    r_outer=2;
end

if nargin <= 2
	points_num = 100;
	bins_num = [8, 10];
end

if nargin <= 1
    nvert = 5;
end

verts = points( round(1 : size(points, 1)/nvert : size(points, 1)), : );

if points_num <= size(points, 1)
    points = points( round(1 : size(points, 1)/points_num : size(points, 1)), : );
end

sc_feat = zeros(nvert, 0);

for i = 1:size(verts,1)
    vertex = verts(i, : ); 
    
    dd = [points(:, 1)-vertex(1), points(:, 2)-vertex(2)];
    % calc radius, and give them radius index
    r = sqrt( dd(:, 1).^2 + dd(:, 2).^2 );
    r = r ./ (mean(r)+eps);
    bin_r = logspace(log10(r_inner), log10(r_outer), bins_num(1)-1);
    bin_r = [0, bin_r];
    cnt_r = zeros(size(r, 1), 1);
    for n = 1 : bins_num(1)
        cnt_r = cnt_r + (r >= bin_r(n));
    end

    % calc angles
    theta = atan2( dd(:, 1), dd(:, 2) ) + pi;
    theta = rem(rem(theta, 2*pi) + 2*pi, 2*pi);
    cnt_theta = 1 + floor( (theta) ./ (2*pi/bins_num(2)) );

    hist = sparse(cnt_r, cnt_theta, 1, bins_num(1), bins_num(2));
    hist = hist(:);
    hist = hist ./ (sum(hist)+eps);
    
    sc_feat(i, 1:length(hist) ) = hist;
end

% vectorize
sc_feat = sc_feat(:);

    

