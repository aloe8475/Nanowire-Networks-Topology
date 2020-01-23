% Make sure you have added the folder 'asn', including its subdiretories
% close all;
% clear;
global glo;
%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;
%% Plot and analysis output flags:
SimulationOptions.takingSnapshots = true; % true \ false
SimulationOptions.compilingMovie  = false; % true \ false 
SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).
%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-2;   % (sec)
SimulationOptions.T  = 240;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  
%% Simulation recording options:
SimulationOptions.ContactMode  = 'farthest';    % '' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
% SimulationOptions.ContactNodes = [randperm(100,glo.contactn)]; % only really required for preSet, other modes will overwrite this
SimulationOptions.ContactNodes = [201,202];
%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
%         Connectivity.filename ='2019-11-18-142407_asn_nw_00100_nj_00876_seed_042_avl_100.00_disp_20.00_gns_05.00_cdisp_80.00.mat';
     Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
%          Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
Connectivity = getConnectivity(Connectivity);
if 0
load('202Mat.mat')
Connectivity.weights= mat202;
Connectivity.EdgeList= edgeList+1;
Connectivity.NumberOfNodes = Connectivity.NumberOfNodes + 2;
Connectivity.NumberOfEdges = Connectivity.NumberOfEdges + 30;
end
%% Choose  contacts:
if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 500; % (um)
end
SimulationOptions = selectContacts(Connectivity, SimulationOptions);
%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);
%% Initialize stimulus:
Stimulus.BiasType              = 'AC';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
switch Stimulus.BiasType
    case 'DC'
        Stimulus.Amplitude = 1.5;  % (Volt)
    case 'AC'
        Stimulus.Frequency       = 0.1; % (Hz)
        Stimulus.Amplitude       = 3;   % (Volt)
    case 'DCandWait'
        Stimulus.OffTime      = 9.9; % SimulationOptions.T/3; % (sec)
        Stimulus.AmplitudeOn  = 1.5;                   % (Volt)
        Stimulus.AmplitudeOff = 0.005;                 % (Volt)
    case 'Ramp'
        Stimulus.AmplitudeMin = 0;    % (Volt)
        Stimulus.AmplitudeMax = 5;    % (Volt)
end
[glo,Stimulus] = getStimuluskevin(glo,Stimulus, SimulationOptions);
%% Get Equations:
% Equations = getEquations(Connectivity,SimulationOptions.ContactNodes);
%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 1*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end
%% Simulate:
if SimulationOptions.takingSnapshots
    [Output, SimulationOptions, snapshots] = simulateNetworknewkevin(glo,Connectivity, Components, Stimulus, SimulationOptions, snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output, SimulationOptions, snapshots] = simulateNetwork(glo,Connectivity, Components, Stimulus, SimulationOptions); % (Ohm)
end
% %% Analysis and plot results:
% if ~SimulationOptions.onlyGraphics
%     plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
% end
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
    switch Stimulus.BiasType
        case 'DC'
            if Stimulus.Signal(1) > 0
                axesLimits.VoltageCbar = [0,Stimulus.Signal(1)]; % (V)
            else
                axesLimits.VoltageCbar = [Stimulus.Signal(1),0]; % (V)
            end
        case {'AC' , 'DCandWait' }
            axesLimits.VoltageCbar = [min(Stimulus.Signal),max(Stimulus.Signal)]; % (V)
    end   
    % Just for fun - extract a specific frame (from the middle):
%    snapshotToFigure(snapshots{floor(length(snapshots)/2)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%     set(gcf, 'visible','on')
    % Compile whole movie:
    if SimulationOptions.compilingMovie 
        fprintf('\nCompiling movie...\n');        
        % Only Windows and MacOS >  10.7
        %v = VideoWriter('networkMovie.mp4','MPEG-4');
        % All platforms
        v = VideoWriter('networkMovie','Motion JPEG AVI');
        v.FrameRate = floor(1/snapshotPeriod/10);       
       v.Quality = 100;
        open(v);
        for i = 1 : length(snapshots)
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