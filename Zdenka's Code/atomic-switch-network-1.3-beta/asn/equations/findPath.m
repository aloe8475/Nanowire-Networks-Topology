function path_ab = findPath(adjMat, v_start, v_finish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds a path from vertex `start` to vertex `finish` in the graph represented 
% by adjMat.
%
% ARGUMENTS: 
% adjMAt -- adjacency matrix of a connected undirected unweighted graph.
% v_start  -- index of the start vertex.
% v_finish -- index of the finish vertex.
%
% OUTPUT:
% path_ab -- a (row) vector of vertex indices specifying a path from start vertex 
% to finish vertex. This path is a sequence of linked nodes/vertices, that never 
% visit a single node more than once. 
%
% REQUIRES:
% none.
%
% Author:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize:
    V       = length(adjMat);
    parent_vertex = zeros(V,1);
    visited_vertex = false(V,1);
    % when visited(i) = true, vertex number i is already in the queue and 
    % its parent node is already known.
    
    % Initialize queue:
    queue = v_start;
    
    % Traverse graph (BFS) until finish vertex is reached:
    while ~visited_vertex(v_finish)
        % Read head of queue:
        currVertex = queue(1);

        % Find all vertices that are adjacent to the current one but have
        % not yet been visited:
        unvisitedAdjacentVertices = adjMat(:,currVertex) & ~visited_vertex;
        
        % Update the parent and flag of these vertices:
        parent_vertex(unvisitedAdjacentVertices)  = currVertex;
        visited_vertex(unvisitedAdjacentVertices) = true;

        % Update queue:
        queue = [queue(2:end); find(unvisitedAdjacentVertices)];
    end

    % Trace back to find a path:
    path_ab = v_finish;
    while path_ab(1) ~= v_start
        path_ab = [parent_vertex(path_ab(1)), path_ab];
    end
end
