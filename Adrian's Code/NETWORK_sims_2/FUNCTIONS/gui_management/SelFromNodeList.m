function [ElSel,NodeSelMat,SelList]=SelFromNodeList(LAdj,NodeList)
LCurr=NodeList.UserData.SelectedNodes;
NList=NodeList.String;
NVal=NodeList.Value;
ElSel=NList(NVal(NVal<=LCurr));
SelList=NList(NVal(NVal>LCurr));
idx=cellfun(@(c) sscanf(c,'Node_%i_%i')',SelList,'UniformOutput',false);
j=cellfun(@(c) c(1),idx);
i=cellfun(@(c) c(2),idx);
NodeSelMat=sparse(i,j,1,LAdj,LAdj);
end