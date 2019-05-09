function varargout=fillNETWORKSETTINGS(ini,inset,handles)
%fill network settings
%shared function between main_sims and callNETWORKSETTINGS

switch ini
    case 0 %modal figure initial
       
        handles.SizeXEdit.String=num2str(inset.SizeX);
        handles.SizeYEdit.String=num2str(inset.SizeY);
        handles.NumberEdit.String=num2str(inset.Number);
        handles.AnglesEdit.String=num2str(inset.Angles);
        handles.LengthEdit.String=num2str(inset.Length);
        handles.ManCheck.Value=inset.ManCheck;
        handles.AngleCheck.Value=inset.AngleCheck;
        handles.DispEdit.String=num2str(inset.Disp);
        handles.NetworkSettings=inset;
        varargout{1}=handles;
      
    case 1 %filling settings form modal
        net.SizeX=str2double(handles.SizeXEdit.String);
        net.SizeY=str2double(handles.SizeYEdit.String);
        net.Number=str2double(handles.NumberEdit.String);
        net.Angles=str2double(handles.AnglesEdit.String);
        net.Length=str2double(handles.LengthEdit.String);
        net.Disp=str2double(handles.DispEdit.String);
        net.ManCheck=handles.ManCheck.Value;
        net.AngleCheck=handles.AngleCheck.Value;
        
        
        varargout{1}=net;
    case 2 %creation of struct in main program
        net.SizeX=20;
        net.SizeY=20;
        net.Number=100;
        net.Angles=88;
        net.Length=6;
        net.Disp=0;
        net.ManCheck=0;
        net.AngleCheck=0;
        
        varargout{1}=net;
    otherwise
        return;
end
end