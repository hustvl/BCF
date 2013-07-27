function [SegmentX, SegmentY, NO]=generate_segments(c, maxvalue, N, nn)

%-------------------------------------------------------
%[SegmentX, SegmentY,NO]=GenSegmentsNew(a,b,maxvalue,nn)
%This function is used to generate all the segments 
%vectors of the input contour
%a and b are the input contour sequence
% maxvalue is the stop condition of DCE, usually 1~1.5
% nn is the sample points' number on each segment, in super's method,n=25

%SegmentX,SegmentY denotes all the coordinates of all the segments of input contour
%NO denotes the number of segments of input contour 
%-------------------------------------------------------

a = c(:,1);
b = c(:,2);

%[x,y,NO]=EvolutionNew2(a,b,ln),
slist=[a,b];

[s, value, delval]=evolution(slist,N,maxvalue,0,0,1);
% [s, value, delval]=evolution(slist,N);

%-------------------------Get the orders of remained Vertexerate_codebook
%of contour.
x=s(:,1);
y=s(:,2);
[t1,t2]=size(s);

ln=t1;

temp=length(a);
j=1;
for i=1:temp;
    if a(i)==x(j)&&b(i)==y(j);
        order(j)=i;
        j=j+1;
        if j>ln;
            break;
        end
       
    end
end
    
%-------------------------Generate every UiUj and UjUi
temp3=1;
scale=ln*(ln-1);
SegmentX=zeros(scale,nn);
SegmentY=zeros(scale,nn);
for i=1:ln;
    for j=i+1:ln;
        temp1=order(i);
        temp2=order(j);
       %if (temp2-temp1)<temp*0.55;
            %----UiUj
        lt=(temp2-temp1)/(nn-1);
     
        T=[0:nn-1]';
        TT=ones(nn,1);
        TTT=round(T*lt+TT*temp1);
        
        
        SegmentX(temp3,:)=a(TTT);
        SegmentY(temp3,:)=b(TTT);
        
        temp3=temp3+1;
        
     % else
        %------- UjUi
        lt=(temp1+temp-temp2)/(nn-1);
       
        T=[0:nn-1]';
        TT=ones(nn,1);
        TTT=round(T*lt+TT*temp2);
        TTTT=find(TTT>temp);
        TTT(TTTT)=TTT(TTTT)-temp;
          
        SegmentX(temp3,:)=a(TTT);
        SegmentY(temp3,:)=b(TTT);
                    
        temp3=temp3+1;
        
      % end
        if temp3>scale;
            break;
        end
    end
end

NO=temp3-1;

SegmentX=SegmentX(1:NO,:);
SegmentY=SegmentY(1:NO,:);

Segments = [SegmentX, SegmentY];