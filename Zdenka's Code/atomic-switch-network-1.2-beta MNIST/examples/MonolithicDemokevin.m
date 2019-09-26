global a;
global plotpsd;
global dutyratio;
global ap;global c;
global b;global M;global t;
global testcurrent;
global contactn;
% Make sure you have added the folder 'asn', including its subdiretories

%% Set the seed for PRNGs for reproducibility:
rng('default');
rng(42)
s = rng;

%% Plot and analysis output flags:
SimulationOptions.takingSnapshots = true; % true \ false
SimulationOptions.compilingMovie  =0; % true \ false 
SimulationOptions.onlyGraphics    = 1; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-2;   % (sec)
SimulationOptions.T  =28*0.1;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Simulation recording options:
SimulationOptions.ContactMode  = 'preSet';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
SimulationOptions.ContactNodes = contactn; % only really required for preSet, other modes will overwrite this

%% Generate Connectivity:

Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
        Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
Connectivity = getConnectivity(Connectivity);

%% Choose  contacts:
if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 2500; % (um)
end
SimulationOptions = selectContacts(Connectivity, SimulationOptions);

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Initialize stimulus:
Stimulus.BiasType              = 'DC';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
switch Stimulus.BiasType
    case 'DC'
        Stimulus.Amplitude = 1.5;  % (Volt)
        dutyratio=length(SimulationOptions.ContactNodes );
    case 'AC'
        Stimulus.Frequency       = 0.1; % (Hz)
        Stimulus.Amplitude       = 3;   % (Volt)
    case 'DCandWait'
        Stimulus.OffTime      = 30; % (sec)
        Stimulus.AmplitudeOn  = 1;                   % (Volt)
        Stimulus.AmplitudeOff = 0.001;                 % (Volt)
    case 'Ramp'
        Stimulus.AmplitudeMin = 0;    % (Volt)
        Stimulus.AmplitudeMax = 5;    % (Volt)
end
Stimulus = getStimulusMNIST(a,Stimulus, SimulationOptions);

%% Get Equations:
Equations = getEquations(Connectivity,SimulationOptions);



%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 1*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end

%% Simulate:
if SimulationOptions.takingSnapshots
%     [Output, SimulationOptions, snapshots] = simulateNetwork(Connectivity,Equations, Components, Stimulus, SimulationOptions, snapshotsIdx); % (Ohm)
[Output, SimulationOptions, snapshots] = simulateNetworknewkevin(Connectivity, Components, Stimulus, SimulationOptions, snapshotsIdx); %
else % this discards the snaphots
    [Output, SimulationOptions, snapshots] =simulateNetworknewkevin(Connectivity,Equations, Components, Stimulus, SimulationOptions); % (Ohm)
end

%% Analysis and plot results:
if ~SimulationOptions.onlyGraphics
%     plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
%     plotResults(Stimulus.Signal./Output.networkCurrent,Output.networkCurrent,Stimulus);
%      plotResults2(Output.networkCurrent,Stimulus);
end

%% Graphics:
if SimulationOptions.takingSnapshots
    % What to plot:
    whatToPlot = struct(...
                        'Nanowires',    true, ...
                        'Contacts',     true, ...
                        'Dissipation',  true, ...
                        'Lambda',       false, ... #can either plot lambda or dissipation or Vdrop
                        'Currents',     true, ...
                        'Voltages',     true, ...
                        'Labels',       false, ...
                        'VDrop',        true, ... 
                        'GraphRep',     true ...
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
        case 'AC'
            axesLimits.VoltageCbar = [min(Stimulus.Signal),max(Stimulus.Signal)]; % (V)
        case'DCandWait'
            axesLimits.VoltageCbar = [0,Stimulus.Signal(1)];  
        case'Ramp'
            axesLimits.VoltageCbar = [0,Stimulus.Signal(1)];
        case 'Triangle'
            axesLimits.VoltageCbar = [0,Stimulus.Signal(1)];
    end
% Just for fun - extract a specific frame (from the middle):
%     figure2=snapshotToFigure(snapshots{floor(length(snapshots)/5)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);    
%    set(gcf, 'visible','on');
%     figure3=snapshotToFigure(snapshots{floor(length(snapshots)*4/10)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);  
%     set(gcf, 'visible','on');
%     figure4=snapshotToFigure(snapshots{floor(length(snapshots)*5/15)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%    set(gcf, 'visible','on');
%     figure5=snapshotToFigure(snapshots{floor(length(snapshots)*9/10)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%   set(gcf, 'visible','on');
%   snapshotToFigure(snapshots{floor(88)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%     set(gcf, 'visible','on');
    
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
        for i = length(snapshots)*1/120 : length(snapshots)*120/120
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

