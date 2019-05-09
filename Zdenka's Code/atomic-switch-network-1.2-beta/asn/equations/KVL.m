function eqMat = KVL(Connectivity,contact,use_parfor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kirchhoff's Voltage Law.
%
% ARGUMENTS: 
% Connectivity -- Structure that contains the adjacency matrix and its
%                 properties. The graph described by the adjacency matrix 
%                 is said to have V vertices and E edges.
% contact -- (2 x 1) row vector vector with the indices of the vertices 
%            between which the external voltage is applied:
%            contact(1)--(V_ext)--contact(2). 
%            The first contact is biased by the external voltage with 
%            respect to the second one. The seocond contact is always 
%            considered as ground (0 V).
% use_parfor -- flag to enable/disable the use of the parfor
%               by default it does not use it.
%
% OUTPUT: 
% eqMat -- |E|-|V|+2 linearly independent KVL equations, represented as a
%          (|E|-|V|+2)x(|E|+1) matrix. This is also termed the circuit 
%          matrix. The KVL cycles are chosen such that there is *exactly one*
%          cycle containing the external source. The equation that
%          originates from this cycle is always placed at the **last** row.
%          The extra column is due to an added tester resistor connected in 
%          series to the external voltage source and to the network.
%
% CONVENTIONS:
% 1. Branches are directed from low to high vertex index.
%    (i) --->---(j) if i < j.
% 2. The orientation of the branch with the tester resistor is from 
%    (i=contact(2)) --->--- (j=contact(1)) regardless whether i < j or 
%    j > i. So, the braanch leaves contact(2) and enters contact(1).
% 3. The quantity that is integrated along each cycle is the potential 
%    gradient (not minus the gradient), so when an edge is traveresed 
%    **against** its orientation, this is represented in the equation as a 
%    **plus** sign, and vice-versa. 
%    This is important only for the last equation.
% 4. The last equation is the only non-homogeneous equation in the
%    set, as it's the only one containing the externally applied voltage. 
%    The externally applied voltage always appears in the RHS without a 
%    sign change. This implies that this cycle must be traversed in a 
%    specific direction (anticlockwise):
%
%    contact(2)-->--through_the_network-->--contact(1)-->--through_the_source-->--contact(2)
% USAGE:
%{
    % Specify one of the test cases
    Connectivity.WhichMatrix   = 'TestCase';
    
    %Load it:
    Connectivity = getConnectivity(Connectivity); 

    eqMatKVL = KVL(Connectivity, [1, 2])
%}
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unless asked to, don't use parallel for:
    if nargin < 3,
        use_parfor = false;
    end
    
    % Initialize:
    V = Connectivity.NumberOfNodes; % Vertices
    E = Connectivity.NumberOfEdges; % Edges
    EV2 = E-V+2;
    eqMat= zeros(EV2,E+1);

    % Find cycles:
    adjMat   = Connectivity.weights;
    edgeList = Connectivity.EdgeList;
    cycles   = cyclebasis(adjMat,'path');
    
    % Add external source cycle:
    cycles{end+1} = [findPath(adjMat, contact(2), contact(1)), V+1];      
    % (contact(2)) ---- [V_ext] ---- (V+1) --- [tester_resistor] ---- (contact(1))
    % Now we are sure that there is exactly one loop containing the
    % source, that it's the last one and that it's traversed in the
    % right direction (as specified in convention #4).
    % When traversed, the added vertex (V+1) is the last one before 
    % closing a cycle.
    % The resulting cycle will be: 
    % [contact(2), ..., through the network, ..., contact(1), V+1, contact(2)]
    
    % Add an edge for the tester resistor listed as edge number E+1: 
    edgeList = [edgeList, [contact(2); V+1]];
    % This edge is between vertices (contact(2)) and (V+1), so if we work 
    % only according to conventions #1,#4, convention #2 is automatically 
    % satisfied.
    
    % Make sure that we have enough cycles and start working:
    if length(cycles) < EV2,
        error(strcat('AtomicSwitchNetworks:', mfilename,':InconsistentCircuitMatrix'), ...
             ['Not enough KVL equations. Have: ' length(cycles) ....
                                        'Need: ' EV2]); 
    else     
        % Convert each cycle to an equation:
        if ~use_parfor,
            for ii = 1 : EV2               
                % Get current cycle and close it:
                thisCycle = cycles{ii};
                thisCycle = [thisCycle thisCycle(1)]; 
                
                % The ends of every edges along the cycle:
                from_vertex = thisCycle(1:end-1);
                to_vertex   = thisCycle(2:end);

                % The orientation of the edges along the cycle:
                cycle_orientation = sign(from_vertex - to_vertex);
                % if from_vertex < to_vertex, then the orientations of the 
                % loop and of the edge are the same: -1.
                % Edge is traversed along the current.
                % Otherwise the edge is traversed in a direction opposite 
                % to its orientation: +1.
                % Edge is traversed against the current.           
            
                % The indices of the edges along the cycle:
                [~, idx] = ismember([from_vertex; to_vertex].' , edgeList.', 'rows');
                edge_idx = idx;
                [~, idx] = ismember([to_vertex; from_vertex].' , edgeList.', 'rows');
                edge_idx = edge_idx + idx;
                
                % The edge from (contact(1)) to (V+1) is the external 
                % source (not an unknown), so it does not have an index.
                % Remove it:
                cycle_orientation(edge_idx==0) = [];
                edge_idx(edge_idx==0) = [];
                
                % Save the equation for the current cycle:
                eqMat(ii,edge_idx) = cycle_orientation;  
            end
        else
            parfor ii = 1 : EV2
                eqMat_v = zeros(1,E+1);
                
                thisCycle = cycles{ii};
                thisCycle = [thisCycle thisCycle(1)]; 
            
                from_vertex = thisCycle(1:end-1);
                to_vertex   = thisCycle(2:end);
            
                cycle_orientation = sign(from_vertex - to_vertex);
            
                [~, idx] = ismember([from_vertex; to_vertex].' , edgeList.', 'rows');
                edge_idx = idx;
                [~, idx] = ismember([to_vertex; from_vertex].' , edgeList.', 'rows');
                edge_idx = edge_idx + idx;
                
                cycle_orientation(edge_idx==0) = [];
                edge_idx(edge_idx==0) = [];
            
                eqMat_v(1, edge_idx) = cycle_orientation;
                eqMat(ii,:) = eqMat_v;
            end
        end
    end
end