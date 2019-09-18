function [Connectivity] = compute_basic_graph_properties(Connectivity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes basic graph properties of the undirected unweighted adjacency 
% matrix. It also plots the average degree distribution as a loglog 
% histogram.
%
% ARGUMENTS: 
% Connectivity -- A structure containing the options. It must contain a 
%                 field .weights, which is the adjacency matrix
% OUTPUT: 
% Connectivity -- the data is returned in the same structure that is passed 
%                 in. Fields added to GraphProperties:
%                 -  .weights - Adjacency matrix.
%                 -  .EdgeList - A 2XE matrix of vertex indices, where 
%                                each column represents an edge. This 
%                                field follows two conventions:
%                                1. edgeList(1,:) < edgeList(2,:)
%                                2. The index of an edge is defined as the 
%                                   index of the column containing the 
%                                   indices of the two vertices it 
%                                   connects.
%
% REQUIRES:
% Brain Connectivity Toolbox 2017-01-15 
% USAGE:
%{

    % Specify a pregenerated nano-wires matrix
    Connectivity.WhichMatrix = 'nanoWires';
    Connectivity.filename    = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    
    %Load it:
    Connectivity = getConnectivity(Connectivity); 

    %Compute stuff
    Connectivity = compute_basic_graph_properties(Connectivity); 

    

%}
% Matlab 2016b
% Authors:
% Paula Sanz-Leon
% Miro Astore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute stuff here
A = Connectivity.weights;

% Node degree
GraphProperties.degrees    = degrees_und(A);
GraphProperties.av_degree = mean(GraphProperties.degrees);
[GraphProperties.min_deg, GraphProperties.min_deg_idx] = min(GraphProperties.degrees);
[GraphProperties.max_deg, GraphProperties.max_deg_idx] = max(GraphProperties.degrees);

% Find node indices that can be used as good candidates for measurment
% points

GraphProperties.node_idx_floor_av_degree = find(GraphProperties.degrees == floor(GraphProperties.av_degree));
GraphProperties.node_idx_ceil_av_degree = find(GraphProperties.degrees == ceil(GraphProperties.av_degree));

% Pairwise shortest paths
GraphProperties.shortest_paths = distance_bin(A);
% The average shortest path length is the characteristic path length of the network.
[lambda,efficiency,ecc,radius,diameter] = charpath(GraphProperties.shortest_paths);

GraphProperties.char_path_length = lambda;
GraphProperties.efficiency_global = efficiency;
GraphProperties.eccentricity_nodal = ecc;
GraphProperties.radius = radius;     % min ecc
GraphProperties.diameter = diameter; % max ecc

Connectivity.GraphProperties = GraphProperties;

end