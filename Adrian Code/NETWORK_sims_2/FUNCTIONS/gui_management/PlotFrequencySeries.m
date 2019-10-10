function PlotFrequencySeries(currAx,NodeList,Sim,type,~,timeInd,varargin)

[ElSel,NodeSelMat]=SelFromNodeList(length(Sim.SelLayout.AdjMat),NodeList);
cla(currAx);
type='-'; %to extract only current from node/electrode extraction functions
cutpoint=50;
if ~isempty(ElSel)
    [T,II,~,~,~]=ExtractCurrentDataFromSim(ElSel,Sim.Data,type,timeInd);
    [F,PSD,beta_legend]=PsdAndBeta(T,II,cutpoint);
    plot(currAx,F,PSD);
    legend(currAx,beta_legend);
    set(currAx,'XScale','log');
    set(currAx,'YScale','log');
end
hold(currAx,'on');
%the price to pay to plot everything
if ~isequal(nnz(NodeSelMat),0)
    [X,Y,~,~,~]=ExtractNodeDataFromSim(NodeSelMat,Sim.Data,type,timeInd);
    [F,PSD,beta_legend]=PsdAndBeta(X,Y,cutpoint);
    plot(currAx,F,PSD);
    legend(currAx,beta_legend);
    set(currAx,'XScale','log');
    set(currAx,'YScale','log');
end

end