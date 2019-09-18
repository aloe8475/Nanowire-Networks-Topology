function fn_pse_single_switch_dc_on_residual(start_idx, end_idx)
    % Quick and dirty parameter space exploration of a single switch.
    % This script runs a simulation for one single atomic switch
    % It saves every single simulation in individual files.
    % Arrays with values of maximum and residual values applied to the switch.
    % The index start_idx and end_idx will access values from this array:
%    dc_amplitude_on  = 2^-9:2^-9:2^-1;       % Volt (step 2^-6: 256 simulations)
%   For a single simulation:
    dc_amplitude_on  = 1:2^-3:1;       % Volt (step 2^-6: 256 simulations)
    
    % Simulations will be run for all the values of residual voltage in
    % this array:
%    dc_bias_residual = 2^-10:2^-10:2^-6;     % Volt (step 2^-10: 16 simulations)
%   For a single simulation:
    dc_bias_residual = 0.0005:2^-9:0.0005;     % Volt (step 2^-10: 16 simulations)

%USAGE:
%{ 
         fn_pse_single_switch_dc_on_residual(1,1) 
         % runs XX simulations for different residual voltage bias. Number of
         % simulations XX depends on length of dc_bias_residual
         
         fn_pse_single_switch_dc_on_residual(1,2) 
         %runs 2 * XX simulations; 2 max applied voltage bias in dc_amplitude_on.
         % Number of simulations XX depends on length of dc_bias_residual

%}

%% Some details of our environment...
%Where is the code
    CodeDir = '..';         % Can be full or relative directory path
    ScriptDir = pwd;        % Get full path to this script
    cd(CodeDir)             % Change to code directory
    FullPathCodeDir = pwd;  % Get full path of CodeDir
    ThisScript = mfilename; % Which script is being run
      
    disp(['Script directory: ' ScriptDir])
    disp(['Code directory: ' FullPathCodeDir])
      
    %When and Where did we start:
    if strcmp(filesep,'/') %on a *nix machine, write machine details log...
       system('uname -a')
       % Add toolbox and its subdirectories to path. Assumes we launch script from
       % examples/ folder
       addpath(genpath('asn/'))
       addpath(genpath('external/'))
    end
    cd(ScriptDir)


if isempty(start_idx)
	pbs_array_idx  = getenv('PBS_ARRAYID');
    pbs_idx        = str2double(pbs_array_idx);
	start_idx = pbs_idx;
    end_idx   = pbs_idx;
end

% The first loops iterates over maximum DC amplitude
for pbs_idx = start_idx:end_idx
    % The second loop iterates over residual voltage
	for idx_int=1:length(dc_bias_residual)
    
    disp(['Script started: ' when()]) 
    disp(['Running: ' ThisScript])

    %% Set the seed for PRNGs for reproducibility:
    rng(42)
    s = rng;

    %% Plot and analysis output flags:
    SimulationOptions.takingSnapshots = true;  % true \ false
    SimulationOptions.compilingMovie  = false; % true \ false 
    SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are gene

    %% Simulation general options:
    SimulationOptions.seed = s;    % save
    SimulationOptions.dt = 2^-10;  % (sec)
    SimulationOptions.T  = 2^6;    % (sec) duration of simulation
    SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
    SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Connectivity
    Connectivity.WhichMatrix    = 'Minimal';
    Connectivity = getConnectivity(Connectivity);

    %% Simulation recording options:
    SimulationOptions.ContactMode  = 'preSet';    
    SimulationOptions.ContactNodes = [1, 2]; 
    
    %% Get equations: (KCL + KVL):
    Equations = getEquations(Connectivity,SimulationOptions.ContactNodes, true);
    disp(['Script finished computing KCL+KVL equations: ' when()])



%% Initialize dynamic components:
    Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
    Components = initializeComponents(Connectivity.NumberOfEdges,Components);

    %% Define stimulus (eg, external input)
    Stimulus.BiasType = 'DCandWait';
    Stimulus.OffTime      = SimulationOptions.T/2;       % (sec)
    Stimulus.AmplitudeOn  = dc_amplitude_on(pbs_idx);    % (Volt) - maximum voltage
%    Stimulus.Amplitude  = dc_amplitude_on(pbs_idx);    % for DC bias
    Stimulus.AmplitudeOff = dc_bias_residual(idx_int);   % (Volt) - residual voltage  
    Stimulus = getStimulus(Stimulus, SimulationOptions);

    %% Initialize snapshot time stamps:r
    if SimulationOptions.takingSnapshots
        snapshotPeriod   = 8*SimulationOptions.dt; % (sec) make it a multiple integer of dt
        snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
        snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
    end

    %% Simulate:
    [Output, SimulationOptions, snapshots] = simulateNetwork(Equations, Components, Stimulus, SimulationOptions, snapshotsIdx); % (Ohm)

    %% Save results to the directory of the invoking script
    % clear unwanted variables
    %close all
%    save([ScriptDir filesep datestr(now, 29) '-' mfilename '_on-' num2str(pbs_idx, '%04d') '_off-' num2str(idx_int, '%04d') '.mat'])
    %% Plot results
    plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
    %% When did we finish:
    disp(['Script ended: ' when()])
end
end
