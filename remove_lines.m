function pnts = remove_lines(pnts)

i_line = [];
for i = 1:length(pnts)
    if is_line(pnts{i})
        i_line = [i_line, i];
    end
end
pnts(i_line) = [];

function f = is_line( ps )

[~, s] = svd( bsxfun(@minus, ps, ps(1,:)) );
s = diag(s);

if s(1)/(s(2)+eps) > 100
    f = 1;
%     figure, plot(ps(:,1), ps(:,2));
else
    f = 0;
end