function [GraphObject,CrossMatX,CrossMatY]=GenerateGraphFromNetwork(Network)
%% Generate graph from network index list of connected nanowires
%% for the moent only topological graph (no weights or directions)

Indexes=Network.IndexInt;
Cx=Network.XInt;
Cy=Network.YInt;
LMat=length(Indexes);

% adjacency matrix

AdjMat=zeros(LMat,LMat);
CrossMatX=zeros(LMat,LMat);
CrossMatY=zeros(LMat,LMat);

for i=1:LMat  
    if ~isempty(Indexes{i})
    AdjMat(i,Indexes{i})=1;
    CrossMatX(i,Indexes{i})=Cx{i};
    CrossMatY(i,Indexes{i})=Cy{i};
    end
end

GraphObject=graph(AdjMat);
end