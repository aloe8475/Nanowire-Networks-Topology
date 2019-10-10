% originally part of cyclebasis.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,E]=spanforest(G)
% determines the spanning forest of G (since may need more than one tree)
% and also the leftover edges.  F is a cell array of adjacency matrices of
% spanning trees, one per disjoint component of G.  E is a cell array of
% adjacency matrices of the complement of the spanning tree wrt the
% corresponding component of G.

F = {};
E = {};

if ~isempty(G),
  spanned = zeros(1,size(G,2)); % vertices which have been spanned          % Logical (binary) vector.

  thisF = zeros(size(G)); % current F                                       % Memory allocation for two adjacency matrices.
  thisE = zeros(size(G)); % current E
  thisv = 1; % vertices of the spanning tree eligible to expand upon        % Always start spanning from vertex number 1.
  spanned(thisv)=1;

  while ~isempty(thisv),
    nextv = [];                                                             % IN THE ORIGINAL CODE LINE (nextv = [];) AND LINE (for i=thisv) WHERE IN REVERSED ORDER. LOGIC MISTAKE.
    for i=thisv,
      % find nodes to extend the spanning tree
      testleaf = find(G(i,:));  
      for j=testleaf,
        if (spanned(j)==1),
          % (i,j) will make tree a loop, so add it to E if it's not in F
          if thisF(i,j)==0,
            thisE(i,j)=1;
            thisE(j,i)=1;
          end;
        else
          % add (i,j) to spanning tree and add j to leaves to test next
          thisF(i,j)=1;
          thisF(j,i)=1;
          spanned(j)=1;
          nextv(end+1)=j;
        end;
      end;
    end;
    thisv = nextv; % update the list of nodes to test next
    if isempty(thisv),
      % no new nodes to test: have a component to add to F & E
      % add it even if it was empty (could have had an isolated self-loop)
      F{end+1} = thisF;
      E{end+1} = thisE;
      thisF = zeros(size(G));
      thisE = zeros(size(G));
      thisv = find(~spanned,1,'first'); % start off a new component
      spanned(thisv)=1;
    end;
  end; % while ~isempty(thisv),
end; % if ~isempty(G),

end % helper function spanforest(...)