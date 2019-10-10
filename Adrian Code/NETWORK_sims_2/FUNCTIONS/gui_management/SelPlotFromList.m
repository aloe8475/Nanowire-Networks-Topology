function PlotSel=SelPlotFromList(handles)
PlotSel=struct();
PList={'1_1','1_2','2_1','2_2'};
for i=1:length(PList)
    Pl=strcat('Plot',PList{i},'List');
    Pp=strcat('Pop',PList{i});
    pops=cellstr(handles.(Pp).String);
    plots=cellstr(handles.(Pl).String);
    PlotSel.Name{i}=pops{handles.(Pp).Value};
    PlotSel.Type{i}=plots{handles.(Pl).Value};
end


end