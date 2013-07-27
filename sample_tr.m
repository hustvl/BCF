function [fea2, lab2] = sample_tr( fea, lab )

labs = unique(lab);
nums = zeros( 1,length(labs) );

for i = 1:length(labs)
    nums(i) = length( find( lab == labs(i) ) );
end

n = min(nums);

sel = [];
for i = 1:length(labs)
    ind = find( lab == labs(i) );
    r = randperm( length(ind) );
    sel = [ sel; ind(r(1:n)) ];
end

fea2 = fea(sel, : );
lab2 = lab(sel);

