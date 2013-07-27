% --------------------------------------------------------------------
% Written by ChengEn Lu, mailto: luchengen@gmail.com, Date: 2008-04-02
% --------------------------------------------------------------------
% INPUT:
% edgelist: input edge list
% nsize: size scale for curvature
% neighbor: effective radius of link area
% curvature_dev: threshold of curvature that could be linked
% --------------------------------------------------------------------
% OUTPUT:
% retlist: edge list after linking, with the same format of
%           edgelist
% --------------------------------------------------------------------
% This function excutes edge relink in an interative rule
% --------------------------------------------------------------------
function linklist =  relinkedge_iter(edgelist, nsize, neighbor, curvature_dev)
iter_prev = 10000;
iter = iter_prev - 1;
linklist = edgelist;
while(iter < iter_prev) & length(linklist) > neighbor
    iter_prev  = iter;
    linklist =  relinkedge(linklist, nsize, neighbor, curvature_dev);
    iter = length(linklist);    
end

return;
