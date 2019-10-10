% Make sure you have added the folder 'asn', including its subdiretories
close all;
clear;

%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;

%% Plot and analysis output flags:
SimulationOptions.takingSnapshots = true; % true \ false
SimulationOptions.compilingMovie  = false; % true \ false 
SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-3;   % (sec)
SimulationOptions.T  = 1;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Simulation recording options:
% SimulationOptions.ContactMode  = 'farthest';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
% SimulationOptions.ContactNodes = [9, 10]; % only really required for preSet, other modes will overwrite this
% For Multiple Electrodes:
SimulationOptions.NumElectrodes = input('How many electrodes would you like to use? \n');
SimulationOptions.ContactMode  = 'preSet';    
SimulationOptions.ContactNodes = [11,12,13,14]; % Set which electrodes to put sources and drains (source, drain, source, drain)
SimulationOptions.isSource = [1,0,1,0]; % Set which electrodes are sources and which are drains. 

%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
          Connectivity.filename='AdriantoZdenka100nw.mat';
%         Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
%         Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
Connectivity = getConnectivity(Connectivity);

%% Choose  contacts:
if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 500; % (um)
end
SimulationOptions = selectContacts(Connectivity, SimulationOptions);


%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponentsAdrian(Connectivity.NumberOfNodes,Components,SimulationOptions);

%% Initialize stimulus:
Stimulus.BiasType              = 'AC';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
switch Stimulus.BiasType
    case 'DC'
        Stimulus.Amplitude = 1.5;  % (Volt)
    case 'AC'
        Stimulus.Frequency       = 1.0; % (Hz)
        Stimulus.Amplitude       = 3;   % (Volt)
    case 'DCandWait'
        Stimulus.OffTime      = 9.9; % SimulationOptions.T/3; % (sec)
        Stimulus.AmplitudeOn  = 1.5;                   % (Volt)
        Stimulus.AmplitudeOff = 0.005;                 % (Volt)
    case 'Ramp'
        Stimulus.AmplitudeMin = 0;    % (Volt)
        Stimulus.AmplitudeMax = 5;    % (Volt)
end

%SimulationOptions.isSource = [1,0,1,0]; % pass this in to getStimulus
Stimulus = getStimulus(Stimulus, SimulationOptions); %NEED TO MODIFY TO GET ALL ELECTRODES

% each column is the voltage of each electrode at each time
Stimulus.Signal(:,2) = zeros(SimulationOptions.NumberOfIterations,1);
Stimulus.Signal(:,3) = Stimulus.Signal(:,1);
Stimulus.Signal(:,4) = zeros(SimulationOptions.NumberOfIterations,1);

%n columns for n electrodes

%% Get Equations:
% Equations = getEquations(Connectivity,SimulationOptions.ContactNodes);
Equations.AdjMat = Connectivity.weights;
Equations.NumberOfEdges = Connectivity.NumberOfEdges;
Equations.NumberOfNodes = Connectivity.NumberOfNodes;

%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 4*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end

%% Simulate:
if SimulationOptions.takingSnapshots
    [Output, SimulationOptions, snapshots] = simulateNetworkAdrian(Equations, Components, Stimulus, SimulationOptions, Connectivity, snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output, SimulationOptions, snapshots] = simulateNetworkAdrian(Equations, Components, Stimulus, SimulationOptions, Connectivity); % (Ohm)
end

% %% Analysis and plot results:
% if ~SimulationOptions.onlyGraphics & SimulationOptions.ContactNodes==2
%     plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
% else
%     plotResultsMulti(Output.networkResistance,Output.sumCurrents,Stimulus);
% end

%% Graphics:
if SimulationOptions.takingSnapshots
    % What to plot:
    whatToPlot = struct(...
                        'Nanowires',    true, ...
                        'Contacts',     true, ...
                        'Dissipation',  true, ...
                        'Currents',     false, ... %need to fix this
                        'Voltages',     true  ...
                        );
    
                    
    % Uniform scales:
    axesLimits.DissipationCbar = [0,5]; % (1pW*10^0 : 1pW*10^5)
    axesLimits.CurrentArrowScaling = 100000000;
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
	
    %% 
    % Just for fun - extract a specific frame (from the middle):
    snapshotToFigure(snapshots{floor(length(snapshots)/2)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
    set(gcf, 'visible','on')

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