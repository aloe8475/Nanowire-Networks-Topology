% This script runs an exemplary simulation of a network with 6 nodes,
% depicted in getConnectivity.
% The dynamics correspond to those of plain resistors.
% The name of the file basically specifies:
%                                          - which network
%                                          - dynamics
%                                          - stimulus

% Approximate runtime: HH:MM , Matlab201Xx, Laptop circa 2015
% Approximate memory:  <?> GB
% Approximate storage: <?> MB 
% From command line:
% matlab -nodesktop -nosplash -r TestNetworkResistiveDCBias

%% 
close all;
clear;

%% Some details of our environment...
%Where is the code
CodeDir = '..';         % Can be full or relative directory path
ScriptDir = pwd;        % Get full path to this script
cd(CodeDir)             % Change to code directory
FullPathCodeDir = pwd;  % Get full path of CodeDir
ThisScript = mfilename; % Which script is being run
  
%When and Where did we start:
disp(['Script started: ' when()]) 
if strcmp(filesep,'/'), %on a *nix machine, write machine details log...
   system('uname -a') 
end

disp(['Running: ' ThisScript])
disp(['Script directory: ' ScriptDir])
disp(['Code directory: ' FullPathCodeDir])

%% Connectivity
clear Connectivity
networkType   = 'TestCase';               
Connectivity.WhichMatrix    = networkType;
Connectivity = getConnectivity(Connectivity);

%% Build the electrical network wrt contact
contact = [1, 2];

%% Get equations: (KCL + KVL):
Equations = getEquations(Connectivity,contact);

%% Set dynamics and initialize components:
Components.ComponentType = 'resistor';              
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Define general parameters for the simulation
% time step
SimulationOptions.dt = 1e-3;                                  % (sec)
% simulation length
SimulationOptions.T  = 1e1;                                   % (sec) 
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Define stimulus (eg, external input)
Stimulus.BiasType = 'DC';
Stimulus.Amplitude = 10;                                      % (Volt)
Stimulus = getStimulus(Stimulus, SimulationOptions);

%% Simulate:
[OutputDynamics, ~, ~] = simulateNetwork(Equations, Components, Stimulus, SimulationOptions); % (Ohm)

%% Save results to the directory of the invoking script
save([ScriptDir filesep mfilename '.mat'])

%% When did we finish:
disp([sprintf('\nScript ended: ') when()])

%% Analyze and plot the results
conductance = 1./OutputDynamics.networkResistance;
disp(conductance(1));