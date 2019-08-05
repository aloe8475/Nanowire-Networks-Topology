function Components = initializeComponentsAdrian(V,Components,SimulationOptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes a structure which holds the characteristics and current state
% of all the electrical elemetns in the network.
%
% ARGUMENTS: 
% E - number of components.
% Components - a structure containing all the options for the components. It
%           must contain a field 'Component Type' which must be one of the
%           following strings:
%           - 'resistor' - passive element
%           - 'memristor' - an element with a charge-dependent resistance
%                           function (memristance) 
%           - 'atomicSwitch' - an element in which switching events are 
%                              driven by voltage.
% 
%
% OUTPUT:
% Components - a struct containing all the properties of the electrical
%              components in the network. {identity, type, voltage, 
%              resistance, onResistance, offResistance} are obligatory 
%              fields, other fields depend on 'componentType'.
%
% REQUIRES:
% none
%
% USAGE:
%{
    Components.ComponentType = 'atomicSwitch'; 
    Components = initializeComponents(Connectivity.NumberOfEdges,Components);
%}
%
% Authors:
% Ido Marcus, Joel Hochstetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %initialises to default values if none are given 
    default.onResistance  = 7.77e-5;
    default.offResistance = 1e-7;
    default.filamentState = 0.0;
    %default.OnOrOff       = 0.0;
    default.setVoltage    = 0.3;%1e-2; %0.3
    default.resetVoltage  = 0.01;%1e-3; %0.01
    default.criticalFlux  = 1e-4;%1e-1; %1e-4
    default.maxFlux       = 0.15;%0.15;
    default.penalty       = 1;%1;
    default.boost         = 1;%1;
    
    fields = fieldnames(default);
    for i = 1:numel(fields)
        if isfield(Components, fields{i}) == 0
            Components.(fields{i}) = default.(fields{i});
        end
    end
        
    
    Components.identity      = [ones(V)]; %; 0];          % ALON 02/08/19 - I DON'T THINK WE NEED THIS 0 for a passive resistor, 1 for an active element
        % If one wants an element to be active with probability p,
        % Components.identity      = [rand(E,1) <= p ; 0];
    Components.type          = Components.ComponentType; % type of active elements ('atomicSwitch' \ 'memristor' \ ...)
    Components.voltage       = zeros(V);             % (Volt)
    Components.resistance    = ones(V)*1e-7;             % (Ohm) (memory allocation)
    Components.onResistance  = [ones(V)*Components.onResistance];   % (Ohm) 1/(12.9 kOhm) = conductance quantum
    Components.offResistance = [ones(V)*Components.offResistance]; %*1e7;   % (Ohm) literature values
    switch Components.ComponentType
        case 'memristor'
            Components.charge         = zeros(V);                                   % (Coulomb)
            % parameters of (... _/-\_/-\ ...) shape:
            Components.lowThreshold   = rand(V)*1e-8;                               % (Coulomb) (1V applied across an OFF-state switch will cause it to open up in about 0.1 sec)
            Components.highThreshold  = (1+rand(V)*1e1) .*Components.lowThreshold;  % (Coulomb)
            Components.period         = (2+rand(V))     .*Components.highThreshold; % (Coulomb) (optional, usually memristance is not a periodic function)
            Components.OnOrOff        = []; % Dummy field only required in atomic switch
        case {'atomicSwitch', 'tunnelSwitch', 'quantCSwitch'}
            % parameters of filament formation\dissociation:
            Components.setVoltage    = sparse(ones(V).*Components.setVoltage); %1e-2;    % (Volt) %% sawtooth: 0.3
            Components.resetVoltage  = sparse(ones(V).*Components.resetVoltage); %1e-3;    % (Volt) %% sawtooth: 0.01
            Components.criticalFlux  = sparse(ones(V).*Components.criticalFlux); %1e-1;  % (Volt*sec)  %% sawtooth: 1e-4
            Components.maxFlux       = sparse(ones(V).*Components.maxFlux); %1.5e-1 % (Volt*sec) %% sawtooth: 0.1
            Components.penalty       = Components.penalty; %10
            Components.boost         = Components.boost; %10
            
            Components.filamentState = sparse(ones(V) .* Components.filamentState);        % (Volt*sec)
            Components.OnOrOff       = sparse(true(V) .* (abs(Components.filamentState) > Components.criticalFlux));

            
        case 'resistor'
            Components.identity      = zeros(V);        % 0 for a passive resistor, 1 for an active element
            Components.OnOrOff       = []; % Dummy field only required in atomic swithc

    end
end