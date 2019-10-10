function result = doesPathExist(adjacencyMatrix, start, finish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines whether there exists a path from "start" to "finish" in the 
% graph represented by "adjacencyMatrix".
%
% ARGUMENTS: 
% adjacencyMatrix - a symmetrix, binary, square matrix representing an
%                   undirected unweighted graph.
% start, finish - indices of vertices in the graph.
%
% OUTPUT:
% result - true if there exists a path from start to finish. false
%          otherwise.
%
% REQUIRES:
% none
%
% Author:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize que and "visited" flag:
visited = false(length(adjacencyMatrix), 1);
que     = start;

% BFS traversal:
while ~isempty(que) && ~visited(finish)
    currentVertex = que(1);
    que = que(2:end);
    visited(currentVertex) = true;
%    size(que)
%    size(find(adjacencyMatrix(:,currentVertex) & ~visited))
    que = [que; find(adjacencyMatrix(:,currentVertex) & ~visited)]; %#ok<AGROW>
end

% Check whether "finish" was reached:
result = visited(finish);