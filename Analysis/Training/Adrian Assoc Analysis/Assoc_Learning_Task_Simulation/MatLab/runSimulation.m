% Make sure you have added the folder 'asn', including its subdiretories
close all;
clear;

tic
%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-3;   % (sec)
SimulationOptions.T  = 1e1;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Simulation recording options:
SimulationOptions.ContactMode     = 'preSet';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
SimulationOptions.electrodes      = [73,30];
SimulationOptions.numOfElectrodes = length(SimulationOptions.electrodes);

%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
        Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
        %Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
Connectivity = getConnectivity(Connectivity);

%% Choose  contacts:
if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 1500; % (um)
end

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Initialize stimulus:
Signals = cell(SimulationOptions.numOfElectrodes,1);

Stimulus1.BiasType       = 'SinglePulse';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
Stimulus1.OnTime         = 0.0; 
Stimulus1.OffTime        = 1.0;
Stimulus1.AmplitudeOn    = 1.4;
Stimulus1.AmplitudeOff   = 0.005;
Signals{1,1} = getStimulus(Stimulus1, SimulationOptions);


% Stimulus2.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
% Signals{2,1} = getStimulus(Stimulus2, SimulationOptions);

Stimulus2.BiasType       = 'SinglePulse';
Stimulus2.OnTime         = 2.0; 
Stimulus2.OffTime        = 3.0;
Stimulus2.AmplitudeOn    = -1.4;
Stimulus2.AmplitudeOff   = -0.005;
Signals{2,1} = getStimulus(Stimulus2, SimulationOptions);





Stimulus3.BiasType       = 'SinglePulse';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
Stimulus3.OnTime         = 0.0; 
Stimulus3.OffTime        = 1.0;
Stimulus3.AmplitudeOn    = 1.4;
Stimulus3.AmplitudeOff   = 0.005;
Signals{3,1} = getStimulus(Stimulus3, SimulationOptions);

Stimulus4.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
Signals{4,1} = getStimulus(Stimulus4, SimulationOptions);

%% Simulate:
fprintf('Running simulation ...')
[Output, SimulationOptions, snapshots] = simulateNetwork(Connectivity, Components, Signals, SimulationOptions); % (Ohm)
fprintf('\n')
% run DataExport.m
toc 
%run ShowMe.m