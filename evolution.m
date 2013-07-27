function [s, value, delval]=evolution(slist,number,maxvalue,keepEndpoints,processUntilConvex,show)
%EVOLUTION(slist, number,<maxvalue>,<keepEndpoints>,<processUntilConvex><show>) 
% discrete curve evolution of slist to n points
%input: slist, number of vertices in resulting list
%       optional: maxvalue: stop criterion,if >0.0 value overrides number
%       optional: keepEndpoints. if set to 1, endpoints are NOT deleted
%       optional: processUntilConvex. if set to 1, evolution continues until
%                 shape is convex
%       optional: show: displayflag
%output: s=simplificated list
%			value= vector containing values of remaining vertices in point-order
%			delval= vector containing values of deleted vertices in order of deletion


% this function is not speed-optimized, it is near a simple brute force
% implementation !

% blocking of vertices is taken into account.

if exist('keepEndpoints')==0 %没有这个输入就如何如何如何
    keepEndpoints=0;
end

if exist('processUntilConvex==0')==0
    processUntilConvex=0;
end

% if exist('show')
%     display=show;
% else
%    display=0;
% end
display = 0;

if (exist('maxvalue')) & maxvalue > 0
   number=3;
else
   maxvalue=inf;
end

if (processUntilConvex>0)
    number=3;
    maxValue = inf;
end

if (number<3) | (length(slist)<=3)
   fprintf('WARNING (evolution): less than 3 vertices\n')
   s=slist;
   return
end

delval=[];
peri=polperimetr([slist;slist(1,:)]);%需要多加个起点作为终点，表示闭合%计算多边形的周长
s=slist;

%initialize value vector (containing the information value of each vertex)
value=zeros(length(s),1);
for i=1:length(s)
	value(i)=relevance(s,i,peri,keepEndpoints); %这个就是relevance的,用的point，就是白老师iccv2009的
end


m=-1;
%mainloop
while 1
   
   if display
   	plotslist(s,1);
      fprintf('value of deleted vertex: %7.5f\n',m);
      pause;
   end
   
   if processUntilConvex & isConvex(s)
       break;                   % => ready
   end
   
   if number >= length(s)
      break;						% => ready
   end
   
   
   [m i]=min(value);				%find min. value
   if m > maxvalue				% => ready 
      break
   end
   
   
   % test blocking
   if m > 0
      bf=blocked(s,i);
      
      % if vertex is blocked, find first non-blocked vertex
      % this procedure is separated from the 'usual min-case'
      % for speed-reasons (sort is needed instead of min)
      if bf
         [rel ind]=sort(value);
         j=2;
      	m=1e16;  %inf.
         
         while j <= length(s)
         	i=ind(j);
         	bf=blocked(s,i);
         	if bf==0
         	   m=rel(j);
         	   break;
            end
            j=j+1;
      	end
      
      	if m>maxvalue   	
         	break;					% => ready.
      	end
      end
   end
   
   %delete vertex
   px=s(i,1);
   py=s(i,2);		% keep coordinates
   s(i,:)=[];
   value(i)=[];   
   delval=[delval; m];
      
   
   %neighbouring vertices
   i0=i-1;
   i1=i;
   if i0==0
   	i0=length(s);
	elseif i1>length(s)
   	i1=1;
   end
   
   
   %neighbouring vertices changed value
   value(i0)=relevance(s,i0,peri,keepEndpoints);
   value(i1)=relevance(s,i1,peri,keepEndpoints);
   
end %while
return

%-------------------------------------------
%RELEVANCE
%returns normalized relevance measure of required index

function v=relevance(s,index,peri,keepEndpoints)

if keepEndpoints
    if index==1 | index == length(s(:,1))
        v=inf;
        return;
    end
end

%vertices
i0=index-1;
i1=index;
i2=index+1;
if i0==0
   i0=length(s);
elseif i2>length(s)
   i2=1;
end

%segments
seg1x=s(i1,1)-s(i0,1);
seg1y=s(i1,2)-s(i0,2);
seg2x=s(i1,1)-s(i2,1);
seg2y=s(i1,2)-s(i2,2);

l1=sqrt(seg1x*seg1x + seg1y*seg1y);
l2=sqrt(seg2x*seg2x + seg2y*seg2y);

%turning angle (0-180)
a=180 - acos(([seg1x seg1y]*[seg2x;seg2y])/l1/l2) * 180 / pi;


%relevance measure
v=a*l1*l2 / (l1 + l2);
v=v/peri;  %normalize
return

%-------------------------------------------------------------
%seglength

function l=seglength(p1x,p1y,p2x,p2y)
dx=p2x-p1x;
dy=p2y-p1y;
l=sqrt(dx*dx + dy*dy);
return

%-------------------------------------------------------------
%blocked

function b=blocked(s,i)

% find neighbouring vertices
i0=i-1;
i1=i+1;
if i0==0
   i0=length(s);
elseif i1>length(s)
   i1=1;
end

%bounding box
minx=min([s(i0,1) s(i,1) s(i1,1)]);
miny=min([s(i0,2) s(i,2) s(i1,2)]);
maxx=max([s(i0,1) s(i,1) s(i1,1)]);
maxy=max([s(i0,2) s(i,2) s(i1,2)]);

% check if any boundary-vertex is inside bounding box
% first create index-set v=s\(i0,i,i1)

if i0<i1
   k=[i1 i i0];
elseif i0>i
   k=[i0 i1 i];
elseif i1<i
   k=[i i0 i1];
end
v=1:length(s);
v(k(1))=[];
v(k(2))=[];
v(k(3))=[];

for k=1:length(v)
   px=s(v(k),1);
   py=s(v(k),2);


	% vertex px,py inside boundary-box ?
	b=0;
	if ((px < minx) | (py < miny) | (px > maxx) | (py > maxy)) == 0
	   %inside, now test triangle
	   a=s(i,:)-s(i0,:);   	%a= i0 to i
	   b=s(i1,:)-s(i,:);		%b= i to i1
	   c=s(i0,:)-s(i1,:);	%c= i1 to i0
	   
	   e0=s(i0,:) - [px py];
	   e1=s(i,:) - [px py];
	   e2=s(i1,:) - [px py];
	   
	   d0=det([a;e0]);
	   d1=det([b;e1]);
	   d2=det([c;e2]);
      
          
	   % INSIDE ?
	   b= ((d0>0) & (d1>0) & (d2>0))  |  ((d0<0) & (d1<0) & (d2<0));
   end
   
   if b
      break
   end
end
return


%-------------------------------------------------------------
% check if shape is convex (or concave)

function con=isConvex(s)

    con = 1;
    if length(s(:,1))<4
        return;
    end
    
    % get direction of first curve
    for i=2:length(s(:,1))-1
        curv=curvatureDirection(s,i);
        if curv ~= 0 break
        end
    end
    
    % check if there's a curve oppositely directed
    for j=i+1:length(s(:,1))-1
        curv1=curvatureDirection(s,j);
        if curv1 ~= 0 & curv1 ~= curv
            con = 0;
            break
        end
    end
    
    return
    

function d=curvatureDirection(s,i)
    a=s(i-1,:);
    b=s(i,:);
    c=s(i+1,:);    
    
    d1=b-a;
    d2=c-b
    
    d=sign(det([d1;d2]));
 return

 
% PERIMETER Perimeter (length) of a polygon.
%	P = PERIMETR(s) Returns perimeter (length) 
%	of a polygon (or a sequence of line segments)
% Calculate length ......................
function  p = polperimetr(s)
p = sum(sqrt(diff(s(:,1)).^2+diff(s(:,2)).^2)); 