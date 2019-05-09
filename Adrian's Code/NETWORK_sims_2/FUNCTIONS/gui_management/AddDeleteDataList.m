function [newcell,newnames,pointvalue]=AddDeleteDataList(action,curr_cell,sel_component)
%% function to add or delte data from a cell array containing those data types
% as well as to create/delete a cell array of str containng the name of those
% data types in the same index order.
switch action
    case 'add'
        if ~isempty(curr_cell)
            pointvalue=length(curr_cell)+1;
            newcell=horzcat(curr_cell,sel_component);
        else
            newcell={sel_component};
            pointvalue=1;
        end
        LNames=length(newcell);
        newnames=cell(LNames,1);
        for i=1:LNames
            newnames{i}=strcat(num2str(i),':',newcell{i}.Name);
        end   
    case {'delete','deleteNet','deleteSim'}
        Delstruct.delete='';
        Delstruct.deleteNet='NetList';
        Delstruct.deleteSim='SimList';
        curr_cell(sel_component)=[];
        newcell=curr_cell;
        if isempty(newcell)
            newnames=Delstruct.(action);
            pointvalue=1;
            return
        end
        LNames=length(newcell);
        newnames=cell(LNames,1);
        for i=1:LNames
            newnames{i}=strcat(num2str(i),':',newcell{i}.Name);
            pointvalue=1;
        end
        
        
        
    otherwise
end
end
