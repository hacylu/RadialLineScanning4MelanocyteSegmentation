% given two pts' coordinate, return the Pts that on the line segment that
% conect these two pts.
% Pts1 is a 1-by-2 vector, where Pts1(1) is the row index, 
% Pts1(2) is the column index
% PtsonLineSegment is a n-by-2 matrix, where n is the number of pts on the
% line segment
% a light version of LgetLineSegmentbyTwoPts, the sampling res is low
function PtsonLineSegment=LgetLineSegmentbyTwoPts_light(Pts1,Pts2,sizeIM)

res=0.05; % the resolution of the pts on line we want

p=polyfit([Pts1(2) Pts2(2)],[Pts1(1) Pts2(1)],1);
% in the case of not a vertical line or horizontal line
if abs(p(1))>20 % in the case of a vertical line
    y=min(Pts1(1),Pts2(1)):res:max(Pts1(1),Pts2(1));
    x=zeros(1,length(y));
    x(:)=Pts1(2);
end
    
if abs(p(1))<0.001% in the case of a horizontal line
    x=min(Pts1(2),Pts2(2)):res:max(Pts1(2),Pts2(2));
    y=zeros(1,length(x));
    y(:)=Pts1(1);
end
    
if abs(p(1))>0.001 && abs(p(1))<20 
    % get the line eqn;
%     p=polyfit([Pts1(2) Pts2(2)],[Pts1(1) Pts2(1)],1);
    % if a line is represented as y=ax+b;then
    a=p(1); b=p(2);
    % get the pts (in image) on the line
    x=min(Pts1(2),Pts2(2)):res:max(Pts1(2),Pts2(2));
    y=a*x+b;
    
end
PtsonLineSegment(:,1)=round(y');
PtsonLineSegment(:,2)=round(x');

%% remove the out of bound values
PtsonLineSegment(PtsonLineSegment(:,1)<min(round(Pts1(1)),round(Pts2(1))),:)=[];
PtsonLineSegment(PtsonLineSegment(:,1)>max(round(Pts1(1)),round(Pts2(1))),:)=[];
PtsonLineSegment(PtsonLineSegment(:,2)>max(round(Pts1(2)),round(Pts2(2))),:)=[];
PtsonLineSegment(PtsonLineSegment(:,2)<min(round(Pts1(2)),round(Pts2(2))),:)=[];

PtsonLineSegment(PtsonLineSegment(:,1)<1,:)=[];
PtsonLineSegment(PtsonLineSegment(:,1)>sizeIM(1),:)=[];
PtsonLineSegment(PtsonLineSegment(:,2)>sizeIM(2),:)=[];
PtsonLineSegment(PtsonLineSegment(:,2)<1,:)=[];

%% remove the redundant point
% codelength=100000000;
% Coded=PtsonLineSegment(:,1)*codelength+PtsonLineSegment(:,2);
% valuni=unique(Coded);
% PtsonLineSegment=[];
% PtsonLineSegment(:,1)=valuni/codelength;
% PtsonLineSegment(:,2)=mod(valuni,codelength);
% PtsonLineSegment=uint8(PtsonLineSegment);
end
