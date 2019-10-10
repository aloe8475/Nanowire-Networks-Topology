function [Connectivity] = getConnectivity(Connectivity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a structure with an adjacency matrix and its properties.
%
% ARGUMENTS: 
% Connectivity -- A structure containing the options. It must contain a 
%                 field .WhichMatrix, whose value is a string specifying 
%                 the structure to be loaded. Options:
%                 - 'nanoWires' - loads the output of 
%                                 generate_nanowires_network.py.
%                                 Required fields:
%                                   .filename (no default)
%                 - 'randAdjMat' - A random connected graph with a 
%                                  specified number of nodes and average 
%                                  degree.
%                                  Required fields:
%                                    .NumberOfNodes
%                                    .AverageDegree (average number of
%                                                    edges per vertex) 
%                                                   (no default).
%                 - 'NearestNeighbor' - A 1D ring (vertex (i) is
%                                       connected to vertices (i-1) and
%                                       (i+1)).
%                                       Required fields:
%                                         .NumberOfNodes
%                 - 'Random' - Same as randAdjMat, but does not guarentee
%                              connectedness, and does support non-binary
%                              and non-symmetric weights. For now, the rest
%                              of the code does not support non-binary and
%                              non-symmetric weights.
%                              Required fields:
%                                .NumberOfNodes
%                                .WeightType
%                                .Symmetric
%                 - Several test cases, see below.
%                 If not specified, the following default values are used:
%                 - .NumberOFNodes = 42
%                 - .Symmetric = true
%                       Possible values = {true, false}; 
%                 - .WeightType = 'binary' 
%                       Possible values = {'binary', 'weighted'}
%
% OUTPUT: 
% Connectivity -- the data is returned in the same structure that is passed 
%                 in. Fields added for all values of .WhichMatrix:
%                 -  .weights - Adjacency matrix.
%                 -  .EdgeList - A 2XE matrix of vertex indices, where 
%                                each column represents an edge. This 
%                                field follows two conventions:
%                                1. edgeList(1,:) < edgeList(2,:)
%                                2. The index of an edge is defined as the 
%                                   index of the column containing the 
%                                   indices of the two vertices it 
%                                   connects.
%                 -  .NumberOfNodes
%                 -  .NumberOfEdges
%                 -  .speed - conduction speed (for time delays in the 
%                             future) in mm/ms or m/s.
%                 -  .dx  
%                 -  .NodeStr - A cell array containing strings for 
%                               labelling each region in the matrix.
%                 -  .VertexPosition - Euclidean coordinates for centre of 
%                                      regions, in micrometers.
%                 -  .wireDistances - Matrix of pairwise Euclidean
%                                     distance.
%                 -  .delay - Matrix of time delays between regions.
%
%                 For Fields added only when .WhichMatrix='nanoWires',
%                 see below.
%
% REQUIRES:
% spanforest
%
% USAGE:
%{
    % Specify a random matrix with 100 nodes
    Connectivity.WhichMatrix   = 'Random';
    Connectivity.NumberOfNodes = 100;

    % Or specify a pregenerated nano-wires matrix
    Connectivity.WhichMatrix = 'nanoWires';
    Connectivity.filename    = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    
    %Load it:
    Connectivity = getConnectivity(Connectivity); 
%}

% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign default values:
    if ~isfield(Connectivity,'NumberOfNodes')
        Connectivity.NumberOfNodes = 42;
    end
    
    if ~isfield(Connectivity,'Symmetric')
        Connectivity.Symmetric = true;
    end
    
    if ~isfield(Connectivity,'WeightType')
        Connectivity.WeightType = 'binary';
    end

% Add fields that are specific to Connectivity types
    switch Connectivity.WhichMatrix    
    %---------------------------------------------------------------------%  

        case 'nanoWires'
            load(Connectivity.filename, 'length_x', 'length_y', 'number_of_wires', 'adj_matrix', 'wire_distances');
            load(Connectivity.filename, 'xc', 'yc', 'xi', 'yi', 'xa', 'ya', 'xb', 'yb');
            
            Connectivity.GridSize       = [length_x,length_y];
            Connectivity.NumberOfNodes  = number_of_wires;
            Connectivity.weights        = adj_matrix;
            Connectivity.wireDistances  = wire_distances;
            Connectivity.VertexPosition = [xc;yc].';
            Connectivity.WireEnds       = [xa;ya;xb;yb].';
            Connectivity.EdgePosition   = [xi;yi].';

    %---------------------------------------------------------------------%  

        case 'randAdjMat' 
            temp = rand(Connectivity.NumberOfNodes) < ...
                   Connectivity.AverageDegree/Connectivity.NumberOfNodes;
            temp = tril(temp,-1);
            adjMat = temp + temp';

            % No empty lines:
            for i=find(sum(adjMat)==0)
                next = mod(i,Connectivity.NumberOfNodes)+1;
                adjMat(i,next) = 1;
                adjMat(next,i) = 1;
            end

            % One connected component:
            [F,~] = spanforest(adjMat);
            if length(F) == 1
                disp('Graph happened to be connected.');
            else
                disp('Graph happened NOT to be connected. Connecting...');

                F{end+1} = F {1};
                for i = 1 : length(F)-1
                    prevTree = F{i};
                    currTree = F{i+1};

                    prevVertex = find(sum(prevTree),1);
                    currVertex = find(sum(currTree),1);

                    adjMat(prevVertex,currVertex)=1;
                    adjMat(currVertex,prevVertex)=1;
                end

                [F,~] = spanforest(adjMat);
                if length(F)==1
                    disp('Success!');
                else
                    disp('Failure!');
                end
            end

            Connectivity.weights = adjMat;

    %---------------------------------------------------------------------%  

        case 'NearestNeighbour'
            Connectivity.weights = diag(ones(1,Connectivity.NumberOfNodes-1),1);
            Connectivity.weights = Connectivity.weights + diag(ones(1,Connectivity.NumberOfNodes-1),-1);
            Connectivity.weights(1,end) = 1;
            Connectivity.weights(end,1) = 1;

    %---------------------------------------------------------------------%  

        case 'Random'
            if strcmp(Connectivity.WeightType, 'binary')
                Connectivity.weights  = (randi([0, 1], Connectivity.NumberOfNodes, Connectivity.NumberOfNodes)).*(~ eye(Connectivity.NumberOfNodes));
            else
               Connectivity.weights  = (rand(Connectivity.NumberOfNodes)).*(~ eye(Connectivity.NumberOfNodes));
            end

            if Connectivity.Symmetric
                temp = tril(Connectivity.weights,-1);
                Connectivity.weights = temp + temp.';
            end

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         3         2         6  
        % o--\/\/\--o--\/\/\--o--\/\/\--o
        % |         4                   |  
        % | _\/\/\_ o _\/\/\-------------
        % |         5                   |
        % | _\/\/\_ o _\/\/\-------------


        case 'TestCase' 
            Connectivity.weights = [0,0,1,1,1,0;
                                    0,0,1,0,0,1;
                                    1,1,0,0,0,0;
                                    1,0,0,0,0,1;
                                    1,0,0,0,0,1;
                                    0,1,0,1,1,0];
            Connectivity.NumberOfNodes = 6;
    %---------------------------------------------------------------------%  

    % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o
        % 1         2         3         4         5         6

        case 'TestLinear'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = diag(ones(Connectivity.NumberOfNodes-1, 1), 1) + ...
                                   diag(ones(Connectivity.NumberOfNodes-1, 1), -1); 

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o
        % 1         2         3 

        case 'TestSeries'
            Connectivity.NumberOfNodes = 3;
            Connectivity.weights = [0 1 0;
                                    1 0 1;
                                    0 1 0];

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        %               3         4
        %     |--\/\/\--o--\/\/\--o--\/\/\--|
        %  1--o                             o--2
        %     |--\/\/\--o--\/\/\--o--\/\/\--|
        %               5         6

        case 'TestParallel'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = [0 0 1 0 1 0;
                                    0 0 0 1 0 1;
                                    1 0 0 1 0 0;
                                    0 1 1 0 0 0;
                                    1 0 0 0 0 1
                                    0 1 0 0 1 0];
                
    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o--\/\/\--o
        % 1         2         3         4         5         6         1

        case 'TestCircular'
            Connectivity.NumberOfNodes = 6;
            Connectivity.weights = diag(ones(Connectivity.NumberOfNodes-1, 1),  1) + ...
                                   diag(ones(Connectivity.NumberOfNodes-1, 1), -1); 
            Connectivity.weights(1, end)   = 1;
            Connectivity.weights(end, 1)   = 1;

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         2
        % o--\/\/\--o
        % | \       |
        % \  \      /
        % /   -     \
        % \    -    /
        % /       \ \
        % |        \|
        % o--\/\/\--o
        % 3         4

        case 'TestSquare'
            Connectivity.NumberOfNodes = 4;
            temp = zeros(Connectivity.NumberOfNodes);
            temp(1  ,2:end) = 1;
            temp(2:3,  end) = 1;
            Connectivity.weights = temp + temp';

    %---------------------------------------------------------------------%  

        % This network is intended for debugging purposes
        % 1         2         3
        % o--\/\/\--o--\/\/\--o 
        % |         |         |
        % <         <         <
        % >         >         >
        % <         <         <
        % |4        |5        |6
        % o--\/\/\--o--\/\/\--o 
        % |         |         |
        % <         <         <
        % >         >         >
        % <         <         <
        % |         |         |
        % o--\/\/\--o--\/\/\--o 
        % 7         8         9

        case 'TestRectangular'
            Connectivity.NumberOfNodes = 9;
            temp_vec = [1 1 0 1 1 0 1 1];
            temp = diag(ones(Connectivity.NumberOfNodes-3, 1), 3) + diag(temp_vec, 1);
            Connectivity.weights = temp + temp';
            
        case 'Minimal'
            Connectivity.NumberOfNodes = 2;          
            Connectivity.weights = [0 1; 1 0];

    %---------------------------------------------------------------------%  

    end
 
    % Generate fields that are common to all matrices
    if ~isfield(Connectivity,'EdgeList')
        [ii, jj] = find(tril(Connectivity.weights)); 
        Connectivity.EdgeList = [jj ii]'; 
        % (2XE matrix, must follow the conventions specified above)
    end
 
    Connectivity.NumberOfEdges = size(Connectivity.EdgeList, 2);
 
    if ~isfield(Connectivity,'speed')
        Connectivity.speed = 1.0; 
    end
    
    if ~isfield(Connectivity,'dx')
        Connectivity.dx = 1.0; 
    end
    
    if ~isfield(Connectivity,'NodeStr')
    % Generate nodes labels
        for ns = 1:Connectivity.NumberOfNodes
            Connectivity.NodeStr{ns} = num2str(ns);
        end
    end

    if ~isfield(Connectivity,'VertexPosition')
        % Generate position
        Connectivity.VertexPosition = [(1:Connectivity.NumberOfNodes).' zeros(Connectivity.NumberOfNodes,2)];
    end

    if ~isfield(Connectivity,'wireDistances')
        % Generate distance and delay matrix
        wireDistances = [0:Connectivity.NumberOfNodes/2  (Connectivity.NumberOfNodes/2-1):-1:1];
        for n=2:Connectivity.NumberOfNodes 
            wireDistances = [wireDistances ; circshift(wireDistances(n-1,:),[0 1])]; 
        end
        Connectivity.wireDistances = wireDistances.*Connectivity.dx;
    end

    if ~isfield(Connectivity,'delay')
        Connectivity.delay =  Connectivity.wireDistances .* Connectivity.speed;
    end
 
end
