% Make sure you have added the folder 'asn', including its subdiretories
function SelSims= runSimulation(contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage)

close all;
% clearvars -except Nodes testcurrent testconduct contactn recurrent

tic
%% Set the seed for PRNGs for reproducibility:
rng(randseed)
s = rng;

%% Plot and analysis output flags:
SimulationOptions.takingSnapshots = false; % true \ false
SimulationOptions.compilingMovie  = false; % true \ false
SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 1e-2;%1e-3;   % (sec)
SimulationOptions.T  = 2;%1e0;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);

%% Simulation recording options:
if strcmp(biasType,'DC')
    SimulationOptions.ContactMode     = 'farthest';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
    % FOR MULTIPLE ELECTRODES NEED TO DEFINE BELOW AND UNCOMMENT BOTH LINES
else
    SimulationOptions.ContactMode     = 'preSet';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
    SimulationOptions.ContactNodes = contactn; % only really required for preSet, other modes will overwrite this
    SimulationOptions.electrodes = contactn; % only really required for preSet, other modes will overwrite this
end
if strcmp(SimulationOptions.ContactMode,'preSet')
    SimulationOptions.numOfElectrodes = length(SimulationOptions.electrodes);
else
    SimulationOptions.numOfElectrodes = 2;
end
%% Choose  contacts:

if ~strcmp(SimulationOptions.ContactMode,'preSet')
    if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
        SimulationOptions.BiProbeDistance = 1500; % (um)
    end
    SimulationOptions = selectContacts(Connectivity, SimulationOptions);
end

%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
        loadpath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Zdenka Networks\';
        Connectivity.filename = [loadpath 'AdriantoZdenka' num2str(numNanowires) 'nw_simulation1.mat'];
        Connectivity.DataType = 'Zdenka';
        %         Connectivity.filename = 'NetworkData/2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
        %Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
    case 'Lattice'
        Connectivity.sizex = 10;
        Connectivity.sizey = 10;
end
Connectivity = getConnectivityMulti(Connectivity);

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponentsMulti(Connectivity.NumberOfEdges,Components);

%% Initialize stimulus:
Signals = cell(SimulationOptions.numOfElectrodes,1);


Stimulus1.BiasType       = biasType;           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'AlonPulse'
if strcmp(biasType,'TimeDelay')
    Stimulus1.AmplitudeOn  = inputVoltage;
    Stimulus1.AmplitudeOff = 0.005;%1e-3;
    Stimulus1.Period       = 3; %period of the short pulses
    Stimulus1.LongWait     = timeDelay*0.05; %Waiting time between the first set and second set of pulses
    Stimulus1.NumPulse1    = 3; %number of pulses before the long wait
    Stimulus1.NumPulse2    = 1; %number of pulses after the long wait
    
elseif strcmp(biasType,'DCandWait')
    Stimulus1.OnTime         = 0.0;
    Stimulus1.OffTime        = SimulationOptions.T/20; %ten pulses
    Stimulus1.AmplitudeOn    = inputVoltage;
    Stimulus1.AmplitudeOff   = 0.005;
elseif strcmp(biasType,'DC')
    % Stimulus1.OnTime         = 0.0;
    % Stimulus1.OffTime        = SimulationOptions.T/20; %ten pulses
    Stimulus1.AmplitudeOn    = inputVoltage;
    Stimulus1.AmplitudeOff   = 0.005;
end
[Signals{1,1}, Stimulus{1}] = getStimulusMulti(Stimulus1, SimulationOptions);




Stimulus2.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
[Signals{2,1}, Stimulus{2}] = getStimulusMulti(Stimulus2, SimulationOptions);

% Stimulus3.BiasType       = 'DCandWait';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
% Stimulus3.OnTime         = 0.0;
% Stimulus3.OffTime        = 1.0;
% Stimulus3.AmplitudeOn    = 1.4;
% Stimulus3.AmplitudeOff   = 0.005;
% Signals{3,1} = getStimulusMulti(Stimulus3, SimulationOptions);
%
% Stimulus3.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'
% Signals{3,1} = getStimulusMulti(Stimulus3, SimulationOptions);

%% Simulate:
fprintf('Running simulation ...')

%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 1*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end

%% Simulate:

SimulationOptions.electrodes      = SimulationOptions.ContactNodes;
if SimulationOptions.takingSnapshots
    [Output, SimulationOptions, snapshots, SelSims,testcurrent,testconduct] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions, biasType, snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output, SimulationOptions, snapshots, SelSims,testcurrent] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions,biasType); % (Ohm)
end


%Convert Zdenka's structure to Adrian's Structure:
SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions,Connectivity,Components,Stimulus);
SelSims.Settings.SigType = Stimulus{1}.BiasType;
SelSims.NumberOfNodes=Connectivity.NumberOfNodes;
fprintf('\n')

%% Save Simulation
% run DataExport.m
end
%run ShowMe.m