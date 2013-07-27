% --------------------------------------------------------------------
% Written by ChengEn Lu, mailto: luchengen@gmail.com, Date: 2008-03-13
% --------------------------------------------------------------------
% INPUT:
% edgelist:     input edge list
% srcList:      source edge under test
% EndPointList: proper edge list could be linked to srcList
% curvature_dev: threshold of curvature that could be linked
% --------------------------------------------------------------------
% OUTPUT:
% strength: strength of linking, i.e. curvature, the smaller the better
% --------------------------------------------------------------------
% ATTENTION:
% EndPointList should be under the format of [*,*,Y,X]
% -------------------------------------------------------------------- 
function [strength] = linkStrength(edgelist,srcList, EndPointList,nsize)
if nargin < 4 
    nsize = 20;
end
[nr,nc] = size(EndPointList);
srcidx0 = srcList(1,1:2);
mapidx = EndPointList(1:nr,1:2);
srcidx = repmat(srcidx0, nr, 1);

%compute linkStrength
%first rerange the sequence of edge that should be connected
edgesrc = edgelist{srcidx0(1)};
if(srcidx0(2) == 0)
    edgesrc = flipud(edgesrc);
end

for i = 1:nr 
    if(mapidx(i,1) == srcidx0(1))
        strength(i) = 1000;  %avoid linking to the edge itself
    else
        edgedst = edgelist{mapidx(i)};
        if(mapidx(i,2) == 1)
            edgedst = flipud(edgedst);
        end
        %eleminate same points
        if(edgedst(1,:) == edgesrc(length(edgesrc),:));
            edgedst = edgedst(2:length(edgedst),:);
        end
        edge2 = [edgesrc',edgedst']';
        nn = min(length(edgesrc),length(edgedst))-1;
        sz = min(nn,nsize);
        gc = curvature(edge2, length(edgesrc), sz);
        strength(i) = abs(gc);
    end
end


