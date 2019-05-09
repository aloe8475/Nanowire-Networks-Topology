function PlotTimeSeries(currAx,NodeList,Sim,type,status,timeInd,varargin)
%% plot time series of network


[ElSel,NodeSelMat,NodeNames]=SelFromNodeList(length(Sim.SelLayout.AdjMat),NodeList);
currT=Sim.Data.time(timeInd);NL=length(findobj(currAx,'Type','line'));
if ~isequal(status,'redraw') || ~isequal(length(ElSel),NL-1) || ~isequal(ElSel',currAx.Legend.String)
    cla(currAx);
    if ~isempty(ElSel)
        [T,II,~,~,yxscale]=ExtractCurrentDataFromSim(ElSel,Sim.Data,type,timeInd);
        plot(currAx,T,II);
        hold(currAx,'on');
        
        set(currAx,'XScale',yxscale.x);
        set(currAx,'YScale',yxscale.y);
    end
    hold(currAx,'on');
    %the price to pay to plot everything
    if ~isequal(nnz(NodeSelMat),0)
        [X,Y,~,~,yxscale]=ExtractNodeDataFromSim(NodeSelMat,Sim.Data,type,timeInd);
        plot(currAx,X,Y);
        hold(currAx,'on');
        
        set(currAx,'XScale',yxscale.x);
        set(currAx,'YScale',yxscale.y);
    end
    
    xl=currT;
    Ylims=currAx.YLim;
    line(currAx,[xl xl],[Ylims(1) Ylims(2)],'Color','r','LineWidth',1.5,'Tag','CurrT');
    %legend
    if length([ElSel ; NodeNames])<12
        %there is no point in plotting a big legend
        lg=findobj(currAx,'Type','line','LineWidth',0.5);
        legend(flipud(lg),[ElSel;NodeNames],'Location','northwest','FontSize',7,...
            'Interpreter','none');
    else
        legend(currAx,'off');
    end
    
else
    h=findobj(currAx,'Tag','CurrT');
    h.XData=[currT currT];
end



hold off;
end
