function Equations = getEquations(Connectivity, SimulationOptions, use_parfor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the abstract (=wihtout actual conductance values) matrix of 
% coefficients, to be used by runSimulation.
%
% ARGUMENTS: 
% Connectivity -- Structure that contains the adjacency matrix and its
%                 properties. The graph described by the adjacency matrix 
%                 is said to have V vertices and E edges.
% contact -- (2 x 1) row vector with the indices of the vertices between
%            which the external voltage is applied:
%            contact(1)--(V_ext)--contact(2). 
%            The first contact is biased by the external voltage with 
%            respect to the second one. The seocond contact is always 
%            considered as ground (0 V).
% use_parfor -- flag to enable/disable the use of the parfor in KVL
%               by default it does not use it.
%
% OUTPUT:
% Equations - a structure containing the following fields:
%             - 'KCLCoeff' - the output of KCL. Matrix of   (V-1) x (E+1);
%             - 'KVLCoeff' - the output of KVL. Matrix of (E-V+2) x (E+1);
%             - 'NumberOfNodes' - number of nodes/vertices in the circuit (=number
%                                 of wires in the network).
%             - 'NumberOfEdges' - number of edges in the circuit (=number
%                                 of wire junctions in the network).
%             - 'ContactNodes' - the indices of the nodes between which the 
%                                source is connected and have been used to 
%                                generate KCL and KVL matrices.
%
%             Together, [KCLCoeff ; KVLCoeff] is an 'abstract' matrix of 
%             coefficients. It's an (E+1)X(E+1) matrix, where the first V-1 
%             rows are from KCL equations and the remaining E-V+2 rows 
%             are from KVL equations, including the last row which is from 
%             the KVL cycle that contains the external voltage source. All 
%             entries are in {-1,0,+1}. The KCL rows should be updated (by 
%             runSimulation) at every iteration, with the current 
%             conductance values, times the sign (hence the 'abstract' in 
%             the name'). The unknowns are:
%             -  First E columns - the voltage drops across components.
%             -  Last column - the voltage drop across the tester resistor.
%
% REQUIRES:
% 	KCL
% 	KVL
%
% USAGE:
%{
    Connectivity.WhichMatrix = 'preGenerated';
    Connectivity.filename = ...;
    Connectivity = getConnectivity(Connectivity);
    contact = [1,2];

    Equations = getEquations(Connectivity,contact)

    eqMat = [Equations.KCLCoeff; Equations.KVLCoeff]
%}
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        use_parfor = false;
    end

    Equations.KCLCoeff      = KCL(Connectivity,SimulationOptions.ContactNodes);
    Equations.KVLCoeff      = KVL(Connectivity,SimulationOptions.ContactNodes);
    Equations.NumberOfNodes = Connectivity.NumberOfNodes;
    Equations.NumberOfEdges = Connectivity.NumberOfEdges;
    Equations.ContactNodes  = SimulationOptions.ContactNodes; 
    
end