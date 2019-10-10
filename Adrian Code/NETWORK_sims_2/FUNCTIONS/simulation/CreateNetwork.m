function NewNet=CreateNetwork(InfoText,NetworkSettings)

%% CREATE NANOWIRES
NewNet=PlaceNanowires(InfoText,NetworkSettings);

%% FIND INTERSECTIONS 
NewNet=FindIntersections(NewNet,NetworkSettings);

end