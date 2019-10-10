function eqMat = KCL(Connectivity, contact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kirchhoff's Current Law.
%
% ARGUMENTS: 
% Connectivity -- Structure that contains the adjacency matrix and its
%                 properties. The graph described by the adjacency matrix 
%                 is said to have V vertices and E edges.
% contact -- (2 x 1) row vector vector with the indices of the vertices 
%            between which the external voltage is applied:
%            contact(1)--(V_ext)--contact(2). 
%            The first contact is biased by the external voltage with 
%            respect to the second one. The second contact is always 
%            considered as ground (0 V). If the bias is positive, the measured 
%            current will be positive. If the bias is negative the measured 
%            current will be negative.
%
% OUTPUT: 
% eqMat -- |V|-1 "abstract" (= without inverse resistance values), linearly 
%          independent KCL equations, represented in a (|V|-1)x(|E|+1) 
%          matrix. The extra column is due to an added tester resistor 
%          connected in series to the external voltage source and to the 
%          network.
%          
%
% CONVENTIONS:
% 1. Branches are directed from low to high vertex index.
%    (i) --->---(j) if i < j.
% 2. The orientation of the branch with the tester resistor is from 
%    i=contact(2) to j=contact(1), regardless whether i < j or j > i.
% 3. The quantity that is integrated on every vertex is the current
%    entering it, so an **entering** branch is represented in the equation 
%    as a **plus** sign, and vice-versa. 
%{
    % Specify one of the test cases
    Connectivity.WhichMatrix   = 'TestCase';
    
    %Load it:
    Connectivity = getConnectivity(Connectivity); 

    eqMatKCL = KCL(Connectivity, [1, 2])
%}
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialize:
    V = Connectivity.NumberOfNodes; % Vertices
    E = Connectivity.NumberOfEdges; % Edges
    eqMat = zeros(V, E);
    
    edgeList = Connectivity.EdgeList;
    % edgeList(1,:) < edgeList(2,:) always because we create the edge list 
    % from the lower triangle of adjMat.
    
    % Set -1 for edges leaving source vertices
    lidx = sub2ind(size(eqMat), edgeList(1,:), 1:E);
    eqMat(lidx) = -1; 
    
    % Set +1 for edges entering target vertices
    lidx = sub2ind(size(eqMat), edgeList(2,:), 1:E);
    eqMat(lidx) = +1; 
    
    % Add external source column:
    extSource = zeros(V,1);
    
    % enters contact(1)
    extSource(contact(1)) = +1;  
    % leaves contact(2)
    extSource(contact(2)) = -1;
    % According to convention #2, current on the tester(source) flows from 
    % contact(2) to contact(1), as it compensates for the current in
    % the network which is from the biased contact(1) to the grounded
    % contact(2).
        
    eqMat = [eqMat, extSource];
    
    % The set of V KCL equations is linearly dependent. Namely, if we sum
    % all the equations except one, we recover the missing one. We should
    % therefore toss one equation away.
    eqMat = eqMat(1:end-1,:);
    
end