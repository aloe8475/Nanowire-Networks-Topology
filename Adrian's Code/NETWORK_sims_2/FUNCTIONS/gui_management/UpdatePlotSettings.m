function UpdatePlotSettings(currAx,PlotName,PlotType)
switch PlotName
    case 'FrequencySeries'
        xlabel('F(Hz)');
        ylabel('PSD');
        currAx.Box='on';
        currAx.XScale='log';
        currAx.YScale='log';
    case 'TimeSeries'
        xlabel('t(s)');
        ylabel(PlotType);
        currAx.Box='on';
    case 'ElectrodeSeries'
        currAx.Box='on';
        
        
        
    otherwise
        currAx.XColor='none';
        currAx.YColor='none';
        currAx.Color=[0.35 0.35 0.35];
end