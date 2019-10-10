function y=cyclebasis(G,form)
% CYCLEBASIS - find a basis for the cycle subspace of a graph/network
% usage: y=cyclebasis(G)
%        y=cyclebasis(G,form)
% 
% This can be useful for obtaining the equations for Kirchhoff's second
% law for complicated networks.  The cycles in the basis (fundamental
% cycle) correspond to a set of linearly independent conservation
% equations.
%
% INPUTS: G, form
%  G is the adjacency matrix of the undirected graph/network, with G(i,j)=1
%   if vertices i and j are connected by a single edge. Graph edges have
%   no direction. G can be full or sparse and real or logical.
%   Non-logical matrices are set to G=(G~=0).
%   Non-symmetric inputs are symmetrised by G=(G+G')>0.
%  form (optional) is a string specifying the form of output:
%   'adj' (default) for the adjacency matrices of the cycles;
%   'path' for the lists of vertices which the cycles pass through.
%
% OUTPUTS: y
%  y is a cell array of fundamental cycles (single loops) that form a basis
%   of the cycle subspace.  All cycles and disjoint unions of cycles are a
%   product (using xor) of fundamental cycles.  If form='adj', the
%   adjacency matrices of the fundamental cycles are given, and any cycle
%   satisfies C = xor(...(xor(y{i1},y{i2}),...),y{ik}) for some i1,...,ik.
%   If form='path' then the cycles are specified as paths,
%   e.g. [1,2,4] for the cycle passing through vertices 1,2,4.
%
% ALGORITHM:
%  For each component of the graph, find a spanning tree, then for each
%  missing edge add the edge to the tree and successively remove all leaves
%  and what remains is a fundamental cycle.  So for a connected graph of v
%  vertices and e edges, there are e-v+1 fundamental cycles.
%
% EXAMPLE: cycle subspace basis of complete graph on 3 vertices (with &
%  without self-loops at every vertex) and 2 copies of same:
%    cyclebasis(ones(3,3),'path')          % returns {1,2,[1,2,3],3}
%    cyclebasis(ones(3,3)-eye(3,3),'path') % returns {[1,2,3]}
%    cyclebasis(blkdiag(ones(3,3)-eye(3,3),ones(3,3)-eye(3,3)),'path')
%                                          % returns {[1,2,3],[4,5,6]}

% Author: Ben Petschel 26/6/2009
%
% Change history:
%  26/6/2009  - first release
%  15/8/2016  - bug in spanforest corrected by Ido Marcus
%  13/10/2016 - spanforest moved to a different file

% set default value of form, if necessary
if (nargin<2) || isempty(form),
  form = 'adj';
end;

% determine form of G for output
spout = issparse(G);

% ensure that G is logicals
G = (G~=0);
  
% symmetrize G
if ~isequal(G,G'),
  G = (G+G')>0;
end;


% find the spanning forest of G and the remaining edges
[F,E] = spanforest(G);                                                      % F is a cell array of adjacency matrices of spanning trees, one per disjoint component of G.  E is a cell array of adjacency matrices of the complement of the spanning tree wrt the corresponding component of G.

% count the number of edges in the graphs in E (including loops)
ny = sum(cellfun(@(x)sum(sum(x+x.*eye(size(x))))/2,E));                     % This is the number of cycles that need to be found. First the number of loops in every component is calculated (cellfun does it in turn for every cell in E). Then the results are summed over.
y = cell(1,ny);

% for each edge in E, add it to a tree in F and remove the leaves
k=1;                                                                        % This is a loop counter.
for i=1:numel(F),                                                           % For every component of the graph.
  thisE = E{i};
  [ei,ej] = find(thisE,1,'first');                                          % The first edge not in the spanning tree of this component.
  while ~isempty(ei),
    thisE(ei,ej)=0;                                                         % Once we find the loop for this edge, we don't to see it every again.
    thisE(ej,ei)=0;
    thisF = F{i};
    thisF(ei,ej) = 1;                                                       % The missing edge is added to the spanning tree (and closes a cycle). 
    thisF(ej,ei) = 1;
    y{k} = removeleaves(thisF);                                             % Then, leaves are removed and only the cycle remains. The result is saves in the cell array which is returned by the function.
    k = k+1;
    [ei,ej] = find(thisE,1,'first');                                        % Next missing edge. If there is none, it will return empty and the loop will terminate.
  end; % while ~isempty(ei),
end; % for i=1:numel(F),

if isequal(form,'path'),                                                    % If necessary, convert adjacency matrices of cycles to 'paths'.
  % convert all fundamental cycles to path form
  for k=1:numel(y),
    y{k} = adj2path(y{k});
  end;
elseif spout,                                                               % If necessary, output as sparse matrices.
  for k=1:numel(y),
    y{k} = sparse(y{k});
  end;
end;

end % main function cyclebasis(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H=removeleaves(G)
% removes all leaves of a graph

H = G+G.*eye(size(G)); % count loops twice
i = find(sum(H)==1); % all vertices that are on a single edge
while ~isempty(i),
  H(i,:) = 0;                                                               % Remove 'current layer' of leaves.
  H(:,i) = 0;
  i = find(sum(H)==1);                                                      % Find 'next layer'.
end;

H = (H~=0);                                                                 % Make sure the result is logial (binary).

end % helper function removeleaves(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=adj2path(G)
% converts the adjacency matrix of a cycle to a path
% (warning: no error checking done to ensure G is a single cycle!)
if all(G(:)==0),
  L = [];
else
  L = zeros(1,sum(sum(G+G.*eye(size(G))))/2); % counts number of vertices   % The path's length is the number of edges.
  if any(~ismember(sum(G+G.*eye(size(G))),[0,2])),
    error('G cannot be a cycle because degrees not all equal 0 or 2');
  end;
  [i,j] = find(G,1,'first');                                                % Find where to start the traversal from.
  n = size(G,2);                                                            % n is the number of vertices.
  L(1)=j;                                                                   % j is the current vertex (also in the code loop).
  i = 0;                                                                    % i is a place holder for the previous vertex.
  for k=2:numel(L),                                                         % k is a counter indicating how many edges the path will contain once the current iteration is done.
    % find the next vertex adjacent to j that does not equal previous
    nextj = find(G(j,:)&((1:n)~=i),1,'first');
    i = j;                                                                  % A new vertex was found, so j is the new "previous" vertex.
    j = nextj;                                                              % Update current vertex.
    L(k) = j;                                                               % Update path.
  end;
end; % if all(G(:)==0), ... else ...

end % helper function adj2path(...)