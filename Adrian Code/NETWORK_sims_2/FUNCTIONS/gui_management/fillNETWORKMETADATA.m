function metadata=fillNETWORKMETADATA(varargin)
separator='_';
arg=varargin{1};

if isstruct(arg)
    parnet=arg;
    size=num2str(parnet.NetworkSettings.SizeX);
    NoW=num2str(parnet.NetworkSettings.Number);
    
    GraphNet=parnet.Graph;
    
    
    %% fill metadata
    % more metadata related to the network is to be added in
    % the future
    
    metadata.Tag=[];
    metadata.Type='kirchoff based simulator';
    metadata.MaxNodes=numnodes(GraphNet);
    metadata.Date=datestr(now,'mm/dd-yyyy');
    metadata.Hour=datestr(now,'HH:MM:SS');
    metadata.Name=strcat('Net',separator,...
        'Sx:',size,separator,...
        'NoW',NoW,separator,...
        metadata.Date,separator...
        ,metadata.Hour,separator);
else
    metadata=arg;
end


end