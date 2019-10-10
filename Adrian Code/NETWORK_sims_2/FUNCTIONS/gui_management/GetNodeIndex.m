function NodeIndex=GetNodeIndex(LayOut,xo,yo)

%% Retrieve all edges positions
x=cell2num(LayOut.XInt);
y=cell2num(LayOut.YInt);

%% Create bounding box
Tolerance=0.5;
xup=xo+Tolerance;
xdown=xo-Tolerance;

%% Find best fit for intersection edge
y_range=y(xdown<x & x<xup);
[~,ic]=min(abs(yo-y_range));
yo=y_range(ic);

%% Trace back the indexes of the nanowire with the edge in common
%% select as electrode the nanowire with the less number of edges
%% just a criteria which does not really affect the output of the sim.
IndexNodes=cellfun(@(y) find(yo==y),LayOut.YInt,'UniformOutput',false);
IndexNanowire=find(~cellfun(@isempty,IndexNodes));
[~,idx]=min(LayOut.NEdges(IndexNanowire));
IndexNanowire=IndexNanowire(idx);
IndexNodes=[IndexNodes{IndexNanowire}];


PosX=LayOut.XInt{IndexNanowire}(IndexNodes);
PosY=LayOut.YInt{IndexNanowire}(IndexNodes);

%% Return NodeIndex struct
NodeIndex.IndexNanowire=IndexNanowire;
NodeIndex.IndexNodes=IndexNodes;
NodeIndex.PosX=PosX;
NodeIndex.PosY=PosY;


end