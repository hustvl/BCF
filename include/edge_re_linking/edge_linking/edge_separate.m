% --------------------------------------------------------------------
% Written by ChengEn Lu, mailto: luchengen@gmail.com, Date: 2008-03-13
% --------------------------------------------------------------------
% INPUT:
% edgelist: input edge list
% nsize: scale size for curvature
% curvature_dev: threshold, introduce a threshold at the points with curvature larger than curvature_dev.
% --------------------------------------------------------------------
% OUTPUT:
% retlist: edge list after processing.
% it follows the same format with edgelist
% --------------------------------------------------------------------
function [retlist] = edge_separate(edgelist, nsize, curvature_dev)
if nargin < 2 
    nsize = 30;
    curvature_dev = 0.8; 
elseif nargin < 3
    curvature_dev = 0.8;
end

minlen = nsize*2+1;
cnt = 0;
list_len = length(edgelist);
for i = 1 : list_len
    seq_now = edgelist{i};
    len = length(seq_now);
    %we don't consider the edges that is shorter than nsize*2+1
    if(len < minlen)
        cnt = cnt + 1;
        retlist{cnt} = seq_now;
    else
        % get curvature 
        gc = globalcurv(seq_now, nsize, 1);
        [ma, idx] = max(abs(gc));
        if(ma > curvature_dev)
            cnt = cnt+1;
            retlist{cnt} = seq_now(1:(idx+nsize),:);
            cnt = cnt+1;
            retlist{cnt} = seq_now((idx+nsize+1): end,:);
        else
            cnt = cnt + 1;
            retlist{cnt} = seq_now;
        end
    end
    
end


return;