function result = onOrOff(Snapshot, Connectivity, Contact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines whether the network is in a collective ON or in a collective
% OFF state. The network's collective state is defined to be "ON" iff there 
% exists a path betwee the two contacts which passes only through "ON" 
% switches.
%
% ARGUMENTS: 
% Snapshot - a struct containing the voltage, resistance, etc. of the
%            electrical components in the network, at a particular 
%            timestamp.
%            The only field actually needed is snapshot.OnOrOff
% Connectivity - a struct containing the networks graph-representation,
%                as returned by 'getConnectivity'.
%                The only field actually needed is connectivity.weights
% Contact - indices of the nanowires (vertices) to which the 
%           external voltage is connected. 
%
% OUTPUT:
% result - true if the collective state of the network is ON. false
%          otherwise.
%
% REQUIRES:
% none
%
% Author:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the edges which correspond to OFF switches:
badPairs = Connectivity.EdgeList(~Snapshot.OnOrOff);
    % Reminder: EdgeList is a 2XE matrix of vertex indices, where each 
    % column represents an edge. The index of an edge is defined as the 
    % index of the corresponding column in this list.

% Get the original adjacency matrix:
adjacencyMatrix = Connectivity.weights;

% Remove the edges which correspond to OFF switches:
adjacencyMatrix(sub2ind(size(adjacencyMatrix),badPairs(1,:),badPairs(2,:))) = 0;
adjacencyMatrix(sub2ind(size(adjacencyMatrix),badPairs(2,:),badPairs(1,:))) = 0;

% Check whether in the modified adjacency matrix the two contacts are in 
% the same connected component:
result = doesPathExist(adjacencyMatrix, contact(1), contact(2));