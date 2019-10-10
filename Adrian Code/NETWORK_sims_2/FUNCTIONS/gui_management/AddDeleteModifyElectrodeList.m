function [newcell,newnames,pointvalue]=AddDeleteModifyElectrodeList(action,curr_cell,sel_component,sel_index)
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
            newnames{i}=newcell{i}.Name;
        end
    case 'delete'
        if abs(length(curr_cell)-1)<2
            newcell=[];
            newnames=[];
            pointvalue=[];
            return;
        end
        
        curr_cell(sel_index)=[];
        newcell=curr_cell;
        LNames=length(newcell);
        newnames=cell(LNames,1);
        for i=1:LNames
            newnames{i}=newcell{i}.Name;
            pointvalue=1;
        end
    case 'modify'
        curr_cell{sel_index}=sel_component;
        pointvalue=sel_index;
        newcell=curr_cell;
        LNames=length(newcell);
        newnames=cell(LNames,1);
        for i=1:LNames
            newnames{i}=newcell{i}.Name;
            pointvalue=1;
        end
        
        
        
        
    otherwise
end
end