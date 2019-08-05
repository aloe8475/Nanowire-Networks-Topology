% Make sure you have added the folder 'asn', including its subdiretories
close all;
clear;

%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;
load('AdriantoZdenka100nw.mat','adj_matrix');
for count = 2:length(adj_matrix)
if count ==1
    break
else
    i=count-1;
end 
%% Plot and analysis output flags:
SimulationOptions(i).takingSnapshots = true; % true \ false
SimulationOptions(i).compilingMovie  = false; % true \ false 
SimulationOptions(i).onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions(i).seed = s;    % save
SimulationOptions(i).dt = 1e-3;   % (sec)
SimulationOptions(i).T  = 4;    % (sec) duration of simulation
SimulationOptions(i).TimeVector = (SimulationOptions(i).dt:SimulationOptions(i).dt:SimulationOptions(i).T)';
SimulationOptions(i).NumberOfIterations = length(SimulationOptions(i).TimeVector);  

%% Simulation recording options:
SimulationOptions(i).ContactMode  = 'preSet';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
SimulationOptions(i).ContactNodes = [1, count]; % only really required for preSet, other modes will overwrite this

%% Generate Connectivity:
Connectivity(i).WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity(i).WhichMatrix
    case 'nanoWires'
        Connectivity(i).filename='AdriantoZdenka100nw.mat';
%         Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
        %Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
temp = getConnectivity(Connectivity(i));
clear Connectivity;
Connectivity(i) = temp;
clear temp;
%% Choose  contacts:
if strcmp(SimulationOptions(i).ContactMode, 'specifiedDistance')
    SimulationOptions(i).BiProbeDistance = 500; % (um)
end
SimulationOptions(i) = selectContacts(Connectivity(i), SimulationOptions(i));

%% Initialize dynamic components:
Components(i).ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
temp= initializeComponents(Connectivity(i).NumberOfEdges,Components(i));
clear Components;
Components(i) = temp;
clear temp;

%% Initialize stimulus:
Stimulus(i).BiasType              = 'DC';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
switch Stimulus(i).BiasType
    case 'DC'
        Stimulus(i).Amplitude = 1.5;  % (Volt)
    case 'AC'
        Stimulus.Frequency       = 1.0; % (Hz)
        Stimulus.Amplitude       = 3;   % (Volt)
    case 'DCandWait'
        Stimulus.OffTime      = 1.9; % SimulationOptions.T/3; % (sec)
        Stimulus.AmplitudeOn  = 1.5;                   % (Volt)
        Stimulus.AmplitudeOff = 0.005;                 % (Volt)
    case 'Ramp'
        Stimulus.AmplitudeMin = 0;    % (Volt)
        Stimulus.AmplitudeMax = 5;    % (Volt)
end

temp= getStimulus(Stimulus(i), SimulationOptions(i));
clear Stimulus;
Stimulus(i) = temp;
clear temp;
%% Get Equations:
Equations(i) = getEquations(Connectivity(i),SimulationOptions(i).ContactNodes);

%% Initialize snapshot time stamps:
if SimulationOptions(i).takingSnapshots
    snapshotPeriod   = 4*SimulationOptions(i).dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions(i).dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions(i).NumberOfIterations;
end

%% Simulate:
if SimulationOptions(i).takingSnapshots
    [Output{i}, temp, snapshot] = simulateNetwork(Equations(i), Components(i), Stimulus(i), SimulationOptions(i), snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output{i}, temp, snapshot] = simulateNetwork(Equations(i), Components(i), Stimulus(i), SimulationOptions(i)); % (Ohm)
end
snapshots{i}=snapshot;
clear SimulationOptions
SimulationOptions(i)=temp;
clear temp 
%% Analysis and plot results:
% if ~SimulationOptions(i).onlyGraphics
%     plotResults(Output{i}.networkResistance,Output{i}.networkCurrent,Stimulus(i));
% end

%% Graphics:
if SimulationOptions(i).takingSnapshots
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
    axesLimits.CurrentArrowSacling = 50; %0.25
    switch Stimulus(i).BiasType
        case 'DC'
            if Stimulus(i).Signal(1) > 0
                axesLimits.VoltageCbar = [0,Stimulus(i).Signal(1)]; % (V)
            else
                axesLimits.VoltageCbar = [Stimulus(i).Signal(1),0]; % (V)
                
            end
        case {'AC' , 'DCandWait' }
            axesLimits.VoltageCbar = [min(Stimulus(i).Signal),max(Stimulus(i).Signal)]; % (V)
    end

%     % Just for fun - extract a specific frame (from the middle):
%     
%     snapshotToFigure(snapshot{floor(length(snapshot)/2)},SimulationOptions(i).ContactNodes,Connectivity(i),whatToPlot,axesLimits);
%     set(gcf, 'visible','on')

    % Compile whole movie:
    if SimulationOptions(i).compilingMovie 
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

fprintf(['\n ' num2str(i) '\n']);
clear snapshot
end 