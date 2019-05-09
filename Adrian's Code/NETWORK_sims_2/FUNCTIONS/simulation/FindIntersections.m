function NewNet=FindIntersections(NewNet,NetSettings)
%%Finding common intersection between nanowires and creating
%% the graph structure for the network, added to the table

%first rough approximation with generous bounding boxes
NumNW=NetSettings.Number;
BoxSize=2*NetSettings.Length;
x1=NewNet.x1;
x2=NewNet.x2;
y1=NewNet.y1;
y2=NewNet.y2;
CrossIndex=cell(NumNW,1);

for i=1:NumNW  
    InfBoundX=x1(i)-BoxSize;
    SupBoundX=x2(i)+BoxSize;
    InfBoundY=min([y1(i) y2(i)])-BoxSize;
    SupBoundY=max([y1(i) y2(i)])+BoxSize;
    
    IdCross=find(InfBoundX<x1 & x2<SupBoundX & InfBoundY<min(y1,y2) & ...
        SupBoundY>max(y1,y2));
    IdCross(IdCross==i)=[];
    CrossIndex{i}=IdCross';
    
end

NewNet.CrossIndex=CrossIndex;
%% find actual crossing positions 
NewNet=ParametricCrosses(NewNet);


end