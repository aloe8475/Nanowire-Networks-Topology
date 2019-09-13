function Components = initializeComponents(E,Components)
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
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   Components.identity      = [ones(E,1) ; 0];          % 0 for a passive resistor, 1 for an active element
        % If one wants an element to be active with probability p,
        % Components.identity      = [rand(E,1) <= p ; 0];
    Components.type          = Components.ComponentType; % type of active elements ('atomicSwitch' \ 'memristor' \ ...)
    Components.voltage       = zeros(E+1,1);             % (Volt)
    Components.resistance    = zeros(E+1,1);             % (Ohm) (memory allocation)
    Components.onResistance  = [ones(E,1) ; 1e-5]*1e4;   % (Ohm) conductance quantum
    Components.offResistance = [ones(E,1) ; 1e-8]*1e7;   % (Ohm) literature values
    switch Components.ComponentType
        case 'memristor'
            Components.charge         = zeros(E+1,1);                                   % (Coulomb)
            % parameters of (... _/-\_/-\ ...) shape:
            Components.lowThreshold   = rand(E+1,1)*1e-8;                               % (Coulomb) (1V applied across an OFF-state switch will cause it to open up in about 0.1 sec)
            Components.highThreshold  = (1+rand(E+1,1)*1e1) .*Components.lowThreshold;  % (Coulomb)
            Components.period         = (2+rand(E+1,1))     .*Components.highThreshold; % (Coulomb) (optional, usually memristance is not a periodic function)
            Components.OnOrOff        = []; % Dummy field only required in atomic switch
        case 'atomicSwitch'
            Components.filamentState = zeros(E+1,1);        % (Volt*sec)
            Components.OnOrOff       = false(E+1,1);
            % parameters of filament formation\dissociation:
            Components.setVoltage    = ones(E+1,1)*1e-2; %1e-2;    % (Volt) %% sawtooth: 0.3
            Components.resetVoltage  = ones(E+1,1)*1e-3; %1e-3;    % (Volt) %% sawtooth: 0.01
            Components.criticalFlux  = ones(E+1,1)*1e-1; %1e-1;  % (Volt*sec)  %% sawtooth: 1e-4
            Components.maxFlux       = ones(E+1,1)*0.15; %1.5e-1 % (Volt*sec) %% sawtooth: 0.1
        case 'resistor'
            Components.identity      = zeros(E+1,1);        % 0 for a passive resistor, 1 for an active element
            Components.OnOrOff       = []; % Dummy field only required in atomic swithc

    end
end