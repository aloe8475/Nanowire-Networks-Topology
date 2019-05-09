function [time,out]=getELECTRODES(handles,Settings,Network)





elinfo=Settings.ElectrodesInfo;
numEl=length(elinfo);

%% if previously used electrodes are to be used again
if ~isequal(0,handles.ResElCheck.Value) && ~isequal(handles.SimList.String,'SimList')
    
    IndexNet=handles.NetList.Value;
    IndexSim=handles.SimList.Value;
    SelSim=handles.Networks{IndexNet}.Simulations{IndexSim};
    if ~isequal(cellfun(@(c) c.Name,elinfo,'UniformOutput',false)',SelSim.Electrodes.Name)
        electrodes=table();
        
    else
        electrodes=SelSim.Electrodes;
    end
else
    electrodes=table();
    % functon to retrieve points x,y in plot to calculate functions with data
    uiwait(msgbox('Choose source and drain electrodes for simulation'));
    
end


% suppress stupid warning message of matlab tables
warning off MATLAB:table:RowsAddedExistingVars

%% Get first electrode position

AxHandle=gca;
for i=1:numEl
    if isequal(handles.ResElCheck.Value,0)
        if i>numEl
            break;
        end
        strinfo=elinfo{i}.Name;strinfo=strcat('Place_',strinfo,' electrode.');
        uiwait(msgbox(strinfo));
        [x1,y1]=ginput(1);
        y_lim=get(AxHandle,'YLim');
        x_lim=get(AxHandle,'XLim');
        hl1=line([x1 x1],[y_lim(1) y_lim(2)],'Color','r');
        hl11=line([x_lim(1) x_lim(2)],[y1 y1],'Color','r');
        ParEl=GetNodeIndex(Network.LayOut,x1,y1);
        electrodes.IndexNanowire(i)=ParEl.IndexNanowire;
        electrodes.IndexNodes(i)=ParEl.IndexNodes;
        electrodes.PosX(i)=ParEl.PosX;
        electrodes.PosY(i)=ParEl.PosY;
        cid=cellfun(@(c) find(c==ParEl.IndexNanowire),Network.Domains,...
            'UniformOutput',false);
        electrodes.DomIndex(i)=find(~cellfun(@isempty,cid));
        electrodes.PosIndex(i)=cid{~cellfun(@isempty,cid)};
        electrodes.Name{i}=elinfo{i}.Name;
        hold on
        scatter(electrodes.PosX(i),electrodes.PosY(i),100,'k','s');
    end
    [time,electrodes.Value{i},electrodes.OpenFlag{i}]=SignalGenerator(elinfo{i});
end
%%Get second electrode position


if ~isequal(length(unique(electrodes.DomIndex)),1)
    uiwait(msgbox('Network does not connect through all the selected electrodes'));
    out=[];
    return;
end


%%confirm

button=questdlg('Are you sure?',...
    'Cancel','Yes, go ahead','No, cancel','Yes, go ahead');
if isequal(handles.ResElCheck.Value,0)
    delete(hl1);delete(hl11);
end


switch button
    case 'Yes, go ahead'
        out=electrodes;
    case 'No, cancel'
        out=[];
end

warning on MATLAB:table:RowsAddedExistingVars

end