% --------------------------------------------------------------------
% Written by ChengEn Lu, mailto: luchengen@gmail.com, Date: 2008-03-13
% --------------------------------------------------------------------
% INPUT:
% edgelist: input edge list
% nsize: scale size for curvature
% neighbor: effective radius of link area
% curvature_dev: threshold of curvature that could be linked
% --------------------------------------------------------------------
% OUTPUT:
% retlist: edge list after linking, with the same format of
%           edgelist
% --------------------------------------------------------------------
function [retlist] = relinkedge(edgelist, nsize, neighbor, curvature_dev)
if nargin < 2 
    nsize = 20;
    neighbor = 5;
    curvature_dev = 0.5; 
elseif nargin < 3
    neighbor = 5;
    curvature_dev = 0.5;
elseif nargin < 4
    curvature_dev = 0.5;
end
len = length(edgelist);
%get end points of each edge with edge index
EndPointList = zeros(len*2, 3);
for i = 1:len
    [nr,nc] = size(edgelist{i});
    EndPointList(i*2-1:i*2,1) = i;                  %indicate which edge this end point belongs to
    EndPointList(i*2-1:i*2,2) = [0,1]';             %0 indicates the start point of an edge, while 1 indicates the end point
    EndPointList(i*2-1, 3:4) = edgelist{i}(1,:);    %start point
    EndPointList(i*2, 3:4) = edgelist{i}(nr,:);     %end point
end

%get end point coordinates from end point list
EndPoint = EndPointList(:,3:4);
lenEp = length(EndPoint);
linkmark = zeros(lenEp,lenEp);
%initialize kdtree
[tmp,tmp,TreeRoot]=kdtree(EndPoint,[]);
% Now we should link edges with small curvature variance
for i = 1:lenEp
    % Needn't consider any points that is marked already
    bm = find(linkmark(i,:)>0);
    if(length(bm) > 0)
        continue;
    end
    [close,dist,idx]=kdnn(TreeRoot,EndPoint(i,:),5);
    pt = find(dist < neighbor);
    if(length(pt) > 1)
        idx = idx(1:length(pt));
        ddx = find(idx ~= i);
        idx = idx(ddx);
        % give the edge that is ready for linking
        readyLinkList = EndPointList(idx,:);
        [nr,nc] = size(readyLinkList);
        srcedge = EndPointList(i,:);
        strength = linkStrength(edgelist, srcedge, readyLinkList, nsize);
        [problink,idxsub] = min(strength);
        minflag = 0;
        % estimate whether problink is the min of all curvature 
        for j = 1:nr
            srcedge = readyLinkList(j,:);
            strengthtmp = linkStrength(edgelist, srcedge, readyLinkList, nsize);
            if(min(strengthtmp) < problink)
                minflag = 1;
                break;
            end
        end
        % we will link edges with problink is lower than 'curvature_dev'
        if(problink < curvature_dev && minflag == 0)
            bm = find(linkmark(idx(idxsub),:)>0); 
            if(length(bm > 0))
                continue;
            end
            linkmark(i, idx(idxsub)) = 1;
            linkmark(idx(idxsub), i) = 1;
        end
    end
end

% Now we have got the marked matrix of linking for each edge
% generate new edges

% first find the edges that is larger than 1
summark = sum(linkmark, 2);
sumedge = summark(1:2:lenEp-1)+ summark(2:2:lenEp);

%like a FAT we give the mark when program linking
cnt = 0;
fat = zeros(lenEp/2,1);
hfn = lenEp/2;
for i = 1:hfn
    if(sumedge(i) == 1 && fat(i) == 0)
        %found a new start edge point.
        %now we start to search for the whole edge
        edge = edgelist{i};
        %check which end point the connect point belongs to
        idx = i*2;
        %if the point is not an end point, it should be a start point
        if(summark(idx) ~= 1)
            idx = idx - 1;
            %if it is a start point, flip the edge vector up down.
            edge = flipud(edge);
        end
        fat(i) = 1;         %this edge is taken
        %get edges links to this edge
        [tmp, idx_nxt] = find(linkmark(idx,:) == 1);
        idx2_nxt = ceil(idx_nxt/2);
        edge_nxt = edgelist{idx2_nxt};
        %while this edge is an intermediate edge
        while(sumedge(idx2_nxt) == 2) 
            %again check which end point the connect point belongs to
            %if it is an end point, flip the edge vector up down
            if(mod(idx_nxt,2) == 0)
                edge_nxt = flipud(edge_nxt);
                idx = idx_nxt-1;    %switch from end point to start point
            else
                idx = idx_nxt+1;    %switch from start point to end point
            end
            %check whether this two edges has same points and gaps
            joint_start_pt = edge(end,:);
            joint_end_pt = edge_nxt(1,:);
            if(joint_start_pt == joint_end_pt);
                edge_nxt = edge_nxt(2:end,:);
                %combine this two edges
                edge = [edge',edge_nxt']';
            elseif(sqrt(sum((joint_start_pt-joint_end_pt).^2)) > 2^0.5)
                %fill the gap between this two edge with line segment
                lineseg = getLineSeq(joint_start_pt, joint_end_pt);
                lineseg = lineseg(2:end-1,:);
                %combine edges
                edge = [edge',lineseg',edge_nxt']';
            else
                %combine this two edges
                edge = [edge',edge_nxt']';
            end
            %mark the fat
            fat(idx2_nxt) = 1;
            %find the mapped point
            [tmp, idx_nxt] = find(linkmark(idx,:) == 1);
            %get index
            idx2_nxt = ceil(idx_nxt/2);
            edge_nxt = edgelist{idx2_nxt};
        end
        %get the last edge of this edge chain
        if(mod(idx_nxt,2) == 0)
            edge_nxt = flipud(edge_nxt);
        end
        %check whether this two edges has same points and gaps
        joint_start_pt = edge(end,:);
        joint_end_pt = edge_nxt(1,:);
        if(joint_start_pt == joint_end_pt);
            edge_nxt = edge_nxt(2:length(edge_nxt),:);
            %combine this two edges
            edge = [edge',edge_nxt']';
        elseif(sqrt(sum((joint_start_pt-joint_end_pt).^2)) > 2^0.5)
            %fill the gap between this two edge with line segment
            lineseg = getLineSeq(joint_start_pt, joint_end_pt);
            lineseg = lineseg(2:end-1,:);
            %combine edges
            edge = [edge',lineseg',edge_nxt']';
        else
            %combine this two edges
            edge = [edge',edge_nxt']';
        end
        fat(idx2_nxt) = 1;
        cnt = cnt+1;
        retlist{cnt} = edge;
    elseif(sumedge(i) == 0 && fat(i) == 0)
        cnt = cnt + 1;
        retlist{cnt} = edgelist{i};
        fat(i) = 1;
    end
end

%There may exists some enclosed loops, now we should fix this situation
for i = 1:hfn
    if(fat(i) == 0)
        if(sumedge ~= 2)   %this should not happen now
            fprintf('This should not happen, stupid program!');
        else
            %found a new start edge point.
            %now we start to search for the whole edge
            edge = edgelist{i};
            %start from the end point
            idx = i*2;
            %the last point should be the start point of this edge chain
            stop_point = idx - 1;
            fat(i) = 1;         %this edge is taken
            %get edges links to this edge
            [tmp, idx_nxt] = find(linkmark(idx,:) == 1);
            idx2_nxt = ceil(idx_nxt/2);
            edge_nxt = edgelist{idx2_nxt};
            while(idx_nxt ~= stop_point)
                %again check which end point the connect point belongs to
                %if it is an end point, flip the edge vector up down
                if(mod(idx_nxt,2) == 0)
                    edge_nxt = flipud(edge_nxt);
                    idx = idx_nxt-1;    %switch from end point to start point
                else
                    idx = idx_nxt+1;    %switch from start point to end point
                end
                %check whether these two edges has same points
                joint_start_pt = edge(end,:);
                joint_end_pt = edge_nxt(1,:);
                if(joint_start_pt == joint_end_pt);
                    edge_nxt = edge_nxt(2:length(edge_nxt),:);
                    %combine this two edges
                    edge = [edge',edge_nxt']';
                elseif(sqrt(sum((joint_start_pt-joint_end_pt).^2)) > 2^0.5)
                    %fill the gap between this two edge with line segment
                    lineseg = getLineSeq(joint_start_pt, joint_end_pt);
                    lineseg = lineseg(2:end-1,:);
                    %combine edges
                    edge = [edge',lineseg',edge_nxt']';
                else
                    %combine these two edges
                    edge = [edge',edge_nxt']';
                end
                %mark the fat
                fat(idx2_nxt) = 1;
                %find the mapped point
                [tmp, idx_nxt] = find(linkmark(idx,:) == 1);
                %get index
                idx2_nxt = ceil(idx_nxt/2);
                edge_nxt = edgelist{idx2_nxt};
            end
            fat(idx2_nxt) = 1;
            cnt = cnt+1;
            retlist{cnt} = edge;
        end
    end
end

%free kdtree
kdtree([],[],TreeRoot); 