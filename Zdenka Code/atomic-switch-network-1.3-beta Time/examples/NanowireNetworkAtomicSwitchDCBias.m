% This script runs an exemplary simulation using a pregenerated network 
% with 100 nodes (wires). The dynamics correspond to those of atomic switches.

% Approximate runtime: 00:02:00 , Matlab2015b, Laptop circa 2015/ 
% Processor: Intel i5 4th gen
% Approximate memory:  < 1 GB (gui) & < 350 MB (command line)
% Approximate storage (no video): 13.4 MB 
% From command line:
% matlab -nodesktop -nosplash -r "NanowireNetworkAtomicSwitchDCBias"

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
if strcmp(filesep,'/') %on a *nix machine, write machine details log...
   system('uname -a')
   % Add toolbox and its subdirectories to path. Assumes we launch script from
   % examples/ folder
   addpath(genpath('asn/'))
end

disp(['Script started: ' when()]) 
disp(['Running: ' ThisScript])
disp(['Script directory: ' ScriptDir])
disp(['Code directory: ' FullPathCodeDir])


%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;

%% Plot and analysis output flags:
SimulationOptions.takingSnapshots = true;  % true \ false
SimulationOptions.compilingMovie  = false; % true \ false 
SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are gene

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-3;   % (sec)
SimulationOptions.T  = 100; %6e1;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Connectivity
Connectivity.WhichMatrix    = 'nanoWires';
Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';             
Connectivity = getConnectivity(Connectivity);

%% Simulation recording options:
SimulationOptions.ContactMode  = 'farthest';   
SimulationOptions.ContactNodes = [20, 10];  % wires indices uesd for applying external current. Only really required for preSet, other modes will overwrite this

SimulationOptions = selectContacts(Connectivity, SimulationOptions);  %ZK

%% Get equations: (KCL + KVL):
Equations = getEquations(Connectivity,SimulationOptions.ContactNodes);

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Define stimulus (eg, external input)
Stimulus.BiasType = 'DC';
Stimulus.Amplitude = 2.0;                             % (Volt)
Stimulus = getStimulus(Stimulus, SimulationOptions);

%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 10*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end

%% Simulate:
[Output, SimulationOptions, snapshots] = simulateNetwork(Equations, Components, Stimulus, SimulationOptions, snapshotsIdx); % (Ohm)


%% Analysis and plot results:
if ~SimulationOptions.onlyGraphics
    plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
end

%% Graphics:
if SimulationOptions.takingSnapshots
    % What to plot:
    whatToPlot = struct(...
                        'Nanowires',    true, ...
                        'Contacts',     true, ...
                        'Dissipation',  true, ...
                        'Currents',     true, ...
                        'Voltages',     true  ...
                        );
    
                    
    % Uniform scales:
    axesLimits.DissipationCbar = [0,5]; % (1pW*10^0 : 1pW*10^5)
    axesLimits.CurrentArrowSacling = 0.25;
    if Stimulus.Signal(1) > 0
        axesLimits.VoltageCbar = [0,Stimulus.Signal(1)]; % (V)
    else
        axesLimits.VoltageCbar = [Stimulus.Signal(1),0]; % (V)
    end
    
    % Just for fun - extract a specific frame (from the middle):
    snapshotToFigure(snapshots{floor(length(snapshots)/2)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
    set(gcf, 'visible','on')

    % Compile whole movie:
    if SimulationOptions.compilingMovie 
        fprintf('\nCompiling movie...\n');
        
        videofilename = strcat(datestr(now, 29), '-nanowire-network-asn-demo');
        v = VideoWriter(videofilename,'Motion JPEG AVI');
        v.FrameRate = floor(1/snapshotPeriod/10);
        v.Quality = 100;
        
        open(v);
        for i = 1:length(snapshots)
            progressBar(i,length(snapshots));
            frameFig = snapshotToFigure(snapshots{i},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
            writeVideo(v,getframe(frameFig));
            close(frameFig);
        end
        close(v);
        
        fprintf('\nDone.\n');
    end
end


fprintf('\n');
%% Save results to the directory of the invoking script
% % clear unwanted variables
% close all
% clear v s frameFig 
% save([ScriptDir filesep datestr(now, 29) '-' mfilename '.mat'])
cd(ScriptDir) 

%% When did we finish:
disp(['Script ended: ' when()])
