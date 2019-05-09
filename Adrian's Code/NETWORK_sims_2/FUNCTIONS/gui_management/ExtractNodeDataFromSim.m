function [X,Y,cx,cy,yxscale]=ExtractNodeDataFromSim(NodeMat,SimData,type,timeInd)
X=SimData.time;
X=repmat(X,[1 nnz(NodeMat)]);
cx=X(timeInd,:);

cc=cellfun(@(c) full(c(NodeMat~=0)),SimData.Currents,'UniformOutput',false);
Y=[cc{:}]';
cv=cellfun(@(c) full(c(NodeMat~=0)),SimData.VoltDif,'UniformOutput',false);
V=[cv{:}]';

yxscale.x='linear';yxscale.y='linear';

switch type
    case 'CurrentLog'
        yxscale.y='log';
    case 'Voltage'
        Y=V;
    case 'Conductance'
        cond=abs(Y)./abs(V);
        cond(isnan(cond) | ~isfinite(cond))=NaN;
        Y=cond;
    case 'ConductanceLog'
        cond=abs(Y)./abs(V);
        cond(isnan(cond) | ~isfinite(cond))=NaN;
        Y=cond;
        yxscale.y='log';
        
    case 'Resistance'
        res=abs(V)./abs(Y);
        res(~isfinite(res))=NaN;
        Y=res;
    case 'ResistanceLog'
        res=abs(V)./abs(Y);
        res(~isfinite(res))=NaN;
        Y=res;
        yxscale.y='log';
    case 'JunctionWidth'
        cc=cellfun(@(c) full(c(NodeMat~=0)),SimData.Wmat,'UniformOutput',false);
        Y=[cc{:}]';
    case 'JunctionWidthLog'
        cc=cellfun(@(c) full(c(NodeMat~=0)),SimData.Wmat,'UniformOutput',false);
        Y=[cc{:}]';
        yxscale.y='log';
        
end
cy=Y(timeInd,:);

end