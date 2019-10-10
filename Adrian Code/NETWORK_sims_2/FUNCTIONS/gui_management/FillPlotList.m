function handles=FillPlotList(handlesin)
handles=handlesin;
List=handles.PlotList;
names=List.Names;
types=List.Types;
PList={'1_1','1_2','2_1','2_2'};
default_vals=[1 1 2 5];
for i=1:length(PList)
    Pl=strcat('Plot',PList{i},'List');
    Pp=strcat('Pop',PList{i});
    selname=cellstr(handles.(Pp).String{handles.(Pp).Value});
    handles.(Pl).String=types{strcmp(names,selname)};
    handles.(Pl).Value=default_vals(i);
end


end