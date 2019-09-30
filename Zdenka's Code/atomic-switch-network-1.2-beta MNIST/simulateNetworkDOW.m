function [OutputDynamics, SimulationOptions, snapshots] = simulateNetworkDOW(Equations,Connectivity, Components, Stimulus, SimulationOptions, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Customized version for Joel by Ruomin.
% Using the nodal analysis method to solve for voltage.
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
% Ruomin Zhu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %% Initialize:
    global testcurrent;global c;global t;
    compPtr         = ComponentsPtr(Components);        % using this matlab-style pointer to pass the Components structure by reference
    niterations     = SimulationOptions.NumberOfIterations; 
    contactNodes    = SimulationOptions.ContactNodes;
    E               = Connectivity.NumberOfEdges;
    V               = Connectivity.NumberOfNodes;
    edgeList        = Connectivity.EdgeList.';
    RHS             = zeros(V+length(contactNodes),1); 
        RHSZeros      = zeros(Equations.NumberOfEdges,1); % the first E entries in the RHS vector.
    testcurrent=[];
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
        LHS          = zeros(V+length(contactNodes), V+length(contactNodes));
        LHS(1:V,1:V) = Gmat;
        
         % Get LHS (matrix) and RHS (vector) of equation:
        LHS = [Equations.KCLCoeff ./ compPtr.comp.resistance(:,ones(Equations.NumberOfNodes-1,1)).' ; ...
               Equations.KVLCoeff];
        RHS = [RHSZeros ; Stimulus.Signal(ii)];

        % Solve equation:
           testsol=LHS\RHS;
            testcurrent=[testcurrent,testsol];
%             lhs=ihs1;
        lhs = sparse(LHS);
        rhs = sparse(RHS);
        sol = lhs\rhs;
        if ii==1
            ihs1=lhs;
        end
       

        tempWireV = sol(1:V);
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
                snapshots{kk} = frame;
                kk = kk + 1;
        end
    end
    
    % Store some important fields
    SimulationOptions.SnapshotsIdx = snapshots_idx; % Save these to access the right time from .TimeVector.

    % Calculate network resistance and save:
    
    OutputDynamics.networkCurrent     = networkCurrent;
    OutputDynamics.wireVoltage        = wireVoltage;
    OutputDynamics.junctionVoltage    = junctionVoltage;
    OutputDynamics.junctionResistance = junctionResistance;
    OutputDynamics.junctionFilament   = junctionFilament;
    
    
%     OutputDynamics.testerVoltage     = testerVoltage;
%     OutputDynamics.networkCurrent    = testerVoltage ./ compPtr.comp.resistance(end);
%     OutputDynamics.networkResistance = (Stimulus.Signal ./ testerVoltage - 1).*compPtr.comp.resistance(end);
   % ZK: also for local values:
    %OutputDynamics.lambda = lambda_vals;
    %OutputDynamics.storevoltage = voltage_vals;
    
end