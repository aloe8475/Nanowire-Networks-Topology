function script_cluster_run_dcwait_duration_simulation(nw_order, contact_case, start_idx, end_idx)
% This script runs a simulation using a pregenerated network 
% of approximately XXXX nodes. The dynamics correspond to those of atomic switches.

% Approximate runtime: XX:XX:XX , Matlab2016b, Laptop circa 2016/ 
% Processor: Intel i7 6th gen
% Approximate memory:  < XX GB (gui) & < XX MB (command line)
% Approximate storage (no video): XX.XX MB 
% From command line:
% matlab -nodesktop -nosplash -r "script_cluster_run_dcwait_simulation"

% TODO:
% + read pbs index to change the DC amplitude (first bifurcation parameter)
% + read pbs index to change the DC duration  (second bifurcation parameter)
% + set max simulation time
% + set max DC bias
% + set cases for at least 10 networks
%  
%% 
%close all;
%clear;

% Get index to switch cases if running several cases on the cluster
%pbs_array_idx  = getenv('PBS_ARRAYID');
%pbs_array_idx   = '1'; % for debugging purposes
%pbs_idx         = str2double(pbs_array_idx);
%dc_amplitude_on = 2^-10:2^-7:2^3;     % Volt (step 2^-7: 1024 simulations)
dc_amplitude_on  = 2^-10:2^-6:2^2;     % Volt (step 2^-6: 256 simulations)
dc_bias_residual = 2^-10:2^-10:2^-4;   % Volt (step 2^-12: 64 simulations)
%dc_off_time     = [2, 4, 8, 16, 32];  % T/dc_off_time  duration of the stimulation
%[bias_level, bias_residual]= meshgrid(dc_amplitude_on, dc_bias_residual);

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
    end_idx = pbs_idx;
end
for pbs_idx = start_idx:end_idx
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
SimulationOptions.T  = 2^7;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Connectivity
Connectivity.WhichMatrix    = 'nanoWires';

    switch nw_order % order of magnitude
        case '100'
            Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';             
        case '1000'
            Connectivity.filename = '2017-03-07-143226_asn_nw_00961_nj_01615_seed_074_avl_15.00_disp_20.00_gns_03.00_cdisp_700.00.mat';             

        case '3000'
            Connectivity.filename = '2017-03-07-171712_asn_nw_02921_nj_04413_seed_090_avl_15.00_disp_20.00_gns_03.00_cdisp_700.00.mat';             

        case '10000'
            Connectivity.filename = '2017-03-07-231345_asn_nw_14050_nj_28830_seed_084_avl_15.00_disp_20.00_gns_03.00_cdisp_700.00.mat';             
     
    end

    Connectivity = getConnectivity(Connectivity);

    % Check if we saved graph properties
    try
        load(Connectivity.filename, 'GraphProperties')
    catch
        warning([mfilename ':: GraphProperties structure had not been saved'])
    end


    if ~exist('GraphProperties', 'var')
        %% compute basic graph properties
        Connectivity = compute_basic_graph_properties(Connectivity);
        GraphProperties = Connectivity.GraphProperties;
        save(Connectivity.filename, 'GraphProperties', '-append')
        disp(['\n Script finished computing basic graph properties: ' when()])
    end

    %% Simulation recording options:
    SimulationOptions.ContactMode  = 'preSet';   

    % Choose something based on the network properties
    switch contact_case

            % Extreme case
        case 'max_min'
            SimulationOptions.ContactNodes = [GraphProperties.max_deg_idx, GraphProperties.min_deg_idx]; 
        case 'max_av'
            % Intermediate case
            SimulationOptions.ContactNodes = [GraphProperties.max_deg_idx, GraphProperties.node_idx_ceil_av_degree(1)]; 
        case 'min_av'
            % Intermediate case
            SimulationOptions.ContactNodes = [GraphProperties.min_deg_idx, GraphProperties.node_idx_ceil_av_degree(1)]; 
        case 'av_av' 
            % Average case
            SimulationOptions.ContactNodes = [GraphProperties.node_idx_ceil_av_degree(2), GraphProperties.node_idx_ceil_av_degree(1)]; 
        case 'preset'
            SimulationOptions.ContactNodes = [9, 10]; 
    
    end

    try 
        load(Connectivity.filename, 'Equations')
    catch
        warning([mfilename ':: Equations structure had not been saved']);
    end

    if  exist('Equations', 'var') && isequal(Equations.ContactNodes, SimulationOptions.ContactNodes)
        disp('something')
    else
        %% Get equations: (KCL + KVL):
        Equations = getEquations(Connectivity,SimulationOptions.ContactNodes, true);
        save(Connectivity.filename, 'Equations', '-append')
        disp(['\n Script finished computing KCL+KVL equations: ' when()])
    end


%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Define stimulus (eg, external input)
Stimulus.BiasType = 'DCandWait';
Stimulus.OffTime      = SimulationOptions.T/2;       % (sec)
Stimulus.AmplitudeOn  = dc_amplitude_on(pbs_idx);    % (Volt)
Stimulus.AmplitudeOff = dc_bias_residual(idx_int);   % (Volt)  
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
save([ScriptDir filesep datestr(now, 29) '-' mfilename '_on-' num2str(pbs_idx, '%04d') '_off-' num2str(idx_int, '%04d') '.mat'])

%% When did we finish:
disp(['Script ended: ' when()])
end
end
