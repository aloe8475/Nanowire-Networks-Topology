function [OutputDynamics, SimulationOptions, snapshots, SelSims] = simulateNetwork(Connectivity, Components, Stimulus, SimulationOptions, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates the network and finds the resistance between the two contacts
% as a function of time.
%
% ARGUMENTS:
% Equations - Structure that contains the (abstract) matrix of coefficients
%             (as documented in getEquations) and the number of nodes in
%             the circuit.
% Components - Structure that contains the component properties. Every
%              field is a (E+1)x1 vector. The extra component is the
%              tester resistor connected in series to the voltage and to
%              the network.
% Stimulus - Structure that contains the details of the external stimulus
%            (time axis details, voltage signal).
% SimulationOptions - Structure that contains general simulation details that are indepedent of
%           the other structures (eg, dt and simulation length);
% varargin - if not empty, contains an array of indidces in which a
%            snapshot of the resistances and voltages in the network is
%            requested. This indices are based on the length of the simulation.
% OUTPUT:
% OutputDynamics -- is a struct with the activity of the network
%                    .networkResistance - the resistance of the network (between the two
%                     contacts) as a function of time.
%                    .networkCurrent - the overall current from contact (1) to contact (2) as a
%                     function of time.
% Simulationoptions -- same struct as input, with updated field names
% snapshots - a cell array of structs, holding the resistance and voltage
%             values in the network, at the requested time-stamps.

% REQUIRES:
% updateComponentResistance
% updateComponentState
%
% USAGE:
%{
    Connectivity = getConnectivity(Connectivity);
    contact      = [1,2];
    Equations    = getEquations(Connectivity,contact);
    Components   = initializeComponents(Connectivity.NumberOfEdges,Components)
    Stimulus     = getStimulus(Stimulus);
    
    OutputDynamics = runSimulation(Equations, Components, Stimulus);
%}
%
% Authors:
% Ido Marcus
% Paula Sanz-Leon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize:
compPtr         = ComponentsPtr(Components);        % using this matlab-style pointer to pass the Components structure by reference
niterations     = SimulationOptions.NumberOfIterations;
contactNodes    = SimulationOptions.ContactNodes;
E               = Connectivity.NumberOfEdges;
V               = Connectivity.NumberOfNodes;
edgeList        = Connectivity.EdgeList.';
RHS             = zeros(V+2,1);
%% Output stuff
wireVoltage        = zeros(niterations, V);
networkCurrent     = zeros(niterations, 1);
junctionVoltage    = zeros(niterations, E);
junctionResistance = zeros(niterations, E);
junctionFilament   = zeros(niterations, E);
%% If snapshots are requested, allocate memory for them:
if ~isempty(varargin)
    snapshots           = cell(size(varargin{1}));
    snapshots_idx       = sort(varargin{1});
else
    nsnapshots          = 10;
    snapshots           = cell(nsnapshots,1);
    snapshots_idx       = ceil(logspace(log10(1), log10(niterations), nsnapshots));
end
kk = 1; % Counter

%% Solve equation systems for every time step and update:
for ii = 1 : niterations
    % Show progress:
    progressBar(ii,niterations);
    
    % Update resistance values:
    updateComponentResistance(compPtr);
    
    componentConductance = 1./compPtr.comp.resistance(1:E);
    % Get LHS (matrix) and RHS (vector) of equation:
    Gmat = zeros(V,V);
    % This line can be written in a more efficient vetorized way.
    % Something like:
    % Gmat(edgeList(:,1),edgeList(:,2)) = componentConductance;
    % Gmat(edgeList(:,2),edgeList(:,1)) = componentConductance;
    
    for i = 1:E
        Gmat(edgeList(i,1),edgeList(i,2)) = componentConductance(i);
        Gmat(edgeList(i,2),edgeList(i,1)) = componentConductance(i);
    end
    
    Gmat         = diag(sum(Gmat, 1)) - Gmat;
    LHS          = zeros(V+2, V+2);
    LHS(1:V,1:V) = Gmat;
    
    % Again, should be able to vectorize
    LHS(V+1, contactNodes(1)) = 1;
    LHS(V+2, contactNodes(2)) = 1;
    LHS(contactNodes(1), V+1) = 1;
    LHS(contactNodes(2), V+2) = 1;
    RHS(V+1) = Stimulus.Signal(ii);
    
    % Solve equation:
    lhs = sparse(LHS);
    rhs = sparse(RHS);
    sol = lhs\rhs;
    
    %Wire Voltage
    tempWireV = sol(1:V);
    %Junction Voltage
    compPtr.comp.voltage(1:E) = tempWireV(edgeList(:,1)) - tempWireV(edgeList(:,2));
    
    
    % Update element fields:
    updateComponentState(compPtr, Stimulus.dt);    % ZK: changed to allow retrieval of local values
    %[lambda_vals(ii,:), voltage_vals(ii,:)] = updateComponentState(compPtr, Stimulus.dt);
    
    wireVoltage(ii,:)        = sol(1:V);
    networkCurrent(ii)       = sol(end);
    junctionVoltage(ii,:)    = compPtr.comp.voltage(1:E);
    junctionResistance(ii,:) = compPtr.comp.resistance(1:E);
    junctionFilament(ii,:)   = compPtr.comp.filamentState(1:E);
    
    % Record the activity of the whole network
    if find(snapshots_idx == ii)
        frame.Timestamp  = SimulationOptions.TimeVector(ii);
        frame.Voltage    = compPtr.comp.voltage;
        frame.Resistance = compPtr.comp.resistance;
        frame.OnOrOff    = compPtr.comp.OnOrOff;
        frame.filamentState = compPtr.comp.filamentState;
        frame.Current=frame.Voltage./frame.Resistance;
        %Convert Zdenka code to Adrian Structure
        %Wire Current
        for currWire=1 : Connectivity.NumberOfNodes
            % Find the indices of edges (=intersections) relevant for this
            % vertex (=wire):
            relevantEdges = find(Connectivity.EdgeList(1,:) == currWire | Connectivity.EdgeList(2,:) == currWire);
            
            % Sort them according to physical location (left-to-right, or
            % if the wire is vertical then bottom-up):
            if Connectivity.WireEnds(currWire,1) ~= Connectivity.WireEnds(currWire,3)
                [~,I] = sort(Connectivity.EdgePosition(relevantEdges,1));
            else
                [~,I] = sort(Connectivity.EdgePosition(relevantEdges,2));
            end
            relevantEdges = relevantEdges(I);
            
            % Calculate the current along each section of the wire:
            direction = ((currWire ~= Connectivity.EdgeList(1,relevantEdges))-0.5)*2;
            % Using the convention that currents are defined to flow
            % from wires with lower index to wires with higher index,
            % and that in the field EdgeList the upper row always
            % contains lower indices.
            wireCurrents{currWire} = cumsum(frame.Current(relevantEdges(1:end)).*direction(1:end)');
            % The first element in wireCurrents is the current in the
            % section between relevantEdge(1) and relevantEdge(2). We
            % assume that between relevantEdge(1) and the closest wire
            % end there's no current. Then, acording to a KCL
            % equation, there's also no current from relevantEdge(end)
            % to the other wire end (that's why the end-1 in the
            % expression).
            % The only two exceptions are the two contacts, where the
            % contact point is defined as the wire end closest to
            % relevantEdge(end)
        end
        %Save Wire Currents and Voltages in Adrian's Structure
        frame.WireVoltage{kk}=getAbsoluteVoltage(frame, Connectivity, SimulationOptions.ContactNodes);
        frame.WireCurrents{kk}=wireCurrents;
        %Save Wire Currents and Voltages is normal Structure
        SelSims.Data.WireCurrents{kk}=wireCurrents;
        SelSims.Data.WireVoltages{kk}= frame.WireVoltage{kk};
        snapshots{kk} = frame;
        kk = kk + 1;
    end
end

% Store some important fields
SimulationOptions.SnapshotsIdx = snapshots_idx; % Save these to access the right time from .TimeVector.

% Store some important fields
SimulationOptions.SnapshotsIdx = snapshots_idx; % Save these to access the right time from .TimeVector.

% Calculate network resistance and save:

OutputDynamics.networkCurrent     = networkCurrent;
OutputDynamics.wireVoltage        = wireVoltage;
OutputDynamics.junctionVoltage    = junctionVoltage;
OutputDynamics.junctionResistance = junctionResistance;
OutputDynamics.junctionFilament   = junctionFilament;
end