function [newcell,pointvalue]=AddDeleteDataListNoName(action,curr_cell,sel_component)
%% function to add or delte data from a cell array containing those data types
% as well as to create/delete a cell array of str containng the name of those
% data types in the same index order.

% Code added by Alon 08/04/2019
switch action    
    case 'add'
        if ~isempty(curr_cell)
            pointvalue=length(curr_cell)+1;
            newcell=horzcat(curr_cell,sel_component);
        else
            newcell={sel_component};
            pointvalue=1;
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
        
        
        
    otherwise
end 
end 