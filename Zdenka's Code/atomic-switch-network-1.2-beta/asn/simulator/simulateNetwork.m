function [OutputDynamics, SimulationOptions, snapshots] = simulateNetwork(Equations, Components, Stimulus, SimulationOptions, varargin)
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
    compPtr       = ComponentsPtr(Components);        % using this matlab-style pointer to pass the Components structure by reference
    niterations   = SimulationOptions.NumberOfIterations; 
    testerVoltage = zeros(niterations, 1);                      % memory allocation for the voltage on the tester resistor as function of time
    RHSZeros      = zeros(Equations.NumberOfEdges,1); % the first E entries in the RHS vector.
    
    %% Use sparse matrices:
    Equations.KCLCoeff = sparse(Equations.KCLCoeff);
    Equations.KVLCoeff = sparse(Equations.KVLCoeff);
    RHSZeros           = sparse(RHSZeros);
    
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
        
        % Get LHS (matrix) and RHS (vector) of equation:
        LHS = [Equations.KCLCoeff ./ compPtr.comp.resistance(:,ones(Equations.NumberOfNodes-1,1)).' ; ...
               Equations.KVLCoeff];
        RHS = [RHSZeros ; Stimulus.Signal(ii)];
        
        % Solve equation:
        compPtr.comp.voltage = LHS\RHS;
                
        % Update element fields:
        updateComponentState(compPtr, Stimulus.dt);                           
        
        % Record tester voltage:
        testerVoltage(ii) = compPtr.comp.voltage(end);
        
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
    OutputDynamics.testerVoltage     = testerVoltage;
    OutputDynamics.networkCurrent    = testerVoltage ./ compPtr.comp.resistance(end);
    OutputDynamics.networkResistance = (Stimulus.Signal ./ testerVoltage - 1).*compPtr.comp.resistance(end);
    
end