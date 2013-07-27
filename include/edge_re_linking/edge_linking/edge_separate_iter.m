function [seplist] = edge_separate_iter(edgelist, nsize, curvature_dev);
iter_prev = 0;
iter = iter_prev+1;
while(iter > iter_prev)
    iter_prev  = iter;
    edgelist = edge_separate(edgelist, nsize, curvature_dev);
    iter = length(edgelist);    
end
seplist = edgelist;