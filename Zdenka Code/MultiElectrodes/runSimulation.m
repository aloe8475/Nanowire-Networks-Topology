% Make sure you have added the folder 'asn', including its subdiretories
function SelSims= runSimulation(SimSettings,contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage)

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
SimulationOptions.T  = SimSettings.SimulationDuration;%1e0;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);

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
end

%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    
    case 'nanoWires'
        computer=getenv('computername');
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                loadpath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Zdenka Networks\';
            case '' %if on linux
                loadpath='/headnode2/aloe8475/CODE/Data/Raw/Networks/Zdenka Networks/';
                %     case 'LAPTOP-S1BV3HR7'
                %         currentPath='D:\alon_\Research\PhD\CODE\Analysis';
                %case '' %--- Add other computer paths (e.g. Mike)
        end
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
    
SimulationOptions = selectContacts(Connectivity, SimulationOptions);

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponentsMulti(Connectivity.NumberOfEdges,Components);

%% Initialize stimulus:
Signals = cell(SimulationOptions.numOfElectrodes,1);
% StimulusSource = cell(SimulationOptions.numOfElectrodes,1);

%Source Type
if SimSettings.numSources>1 %if multiple types
    StimulusSource(1).BiasType       = biasType{1};           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'AlonPulse'
    StimulusSource(2).BiasType       = biasType{2};           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'AlonPulse'
else
    StimulusSource.BiasType       = biasType;           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'AlonPulse'
end
%Drain Type
StimulusDrain.BiasType       = 'Drain';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp'

for bias=1:SimSettings.numSources %for each stimulus source, select the bias type
    
    switch biasType
        
        case 'TimeDelay'
            StimulusSource.AmplitudeOn  = inputVoltage;
            StimulusSource.AmplitudeOff = 0.005;%1e-3;
            StimulusSource.Period       = 2; %period of the short pulses
            StimulusSource.LongWait     = timeDelay; %Waiting time between the first set and second set of pulses
            StimulusSource.NumPulse1    = 3; %number of pulses before the long wait
            StimulusSource.NumPulse2    = 1; %number of pulses after the long wait
        case 'TimeDelay1'
            if size(biasType,1) >1
                StimulusSource(bias).AmplitudeOn  = inputVoltage;
                StimulusSource(bias).AmplitudeOff = 0.005;%1e-3;
                StimulusSource(bias).Period       = 2; %period of the short pulses
                StimulusSource(bias).WaitTime     = timeDelay; %Waiting time between the first set and second set of pulses
                StimulusSource(bias).NumPulse1    = 1; %number of pulses before the long wait
                StimulusSource(bias).NumPulse2    = 0; %number of pulses after the long wait
            end
        case 'TimeDelay2'
            if size(biasType,1) >1
                StimulusSource(bias).AmplitudeOn  = inputVoltage;
                StimulusSource(bias).AmplitudeOff = 0.005;%1e-3;
                StimulusSource(bias).Period       = 2; %period of the short pulses
                StimulusSource(bias).StartWait     = timeDelay; %Waiting time between the first set and second set of pulses
                StimulusSource(bias).StartTime     = (StimulusSource(bias).Period) + StimulusSource(bias).StartWait; %start after stimulus 1 + wait time
                StimulusSource(bias).NumPulse1    = 0; %number of pulses before the long wait
                StimulusSource(bias).NumPulse2    = 1; %number of pulses after the long wait
            end
        case 'DCandWait'
            numPulses=10;
            StimulusSource.OnTime         = 0.0;
            StimulusSource.OffTime        = SimulationOptions.T/(numPulses*2); %ten pulses
            StimulusSource.AmplitudeOn    = inputVoltage;
            StimulusSource.AmplitudeOff   = 0.005;
            StimulusSource.Period      = numPulses;
            
        case 'DC'
            % Stimulus1.OnTime         = 0.0;
            % Stimulus1.OffTime        = SimulationOptions.T/20; %ten pulses
            StimulusSource.AmplitudeOn    = inputVoltage;
            StimulusSource.AmplitudeOff   = 0.005;
    end
end
%Source1
[Signals{1,1}, Stimulus{1}] = getStimulusMulti(StimulusSource(1), SimulationOptions);
%Drain1
[Signals{2,1}, Stimulus{2}] = getStimulusMulti(StimulusDrain, SimulationOptions);

if length(StimulusSource)>1%if we want multiple electrodes
    %Source2
    [Signals{3,1}, Stimulus{3}] = getStimulusMulti(StimulusSource(2), SimulationOptions);
    %Drain2
    [Signals{4,1}, Stimulus{4}] = getStimulusMulti(StimulusDrain, SimulationOptions);
end

%% Initialize snapshot time stamps:
if SimulationOptions.takingSnapshots
    snapshotPeriod   = 1*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotStep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotStep:SimulationOptions.NumberOfIterations;
end

%% Simulate:

SimulationOptions.electrodes      = SimulationOptions.ContactNodes;
if SimulationOptions.takingSnapshots
    [Output, SimulationOptions, snapshots, SelSims] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions, biasType, snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output, SimulationOptions, snapshots, SelSims] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions,biasType); % (Ohm)
end


%Convert Zdenka's structure to Adrian's Structure:
SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions,Connectivity,Components,Stimulus,Signals);
if length(StimulusSource)>1
    SelSims.Settings.SigType{1} = Stimulus{1}.BiasType;
    SelSims.Settings.SigType{2} = Stimulus{2}.BiasType;
else
    SelSims.Settings.SigType = Stimulus{1}.BiasType;
end
SelSims.NumberOfNodes=Connectivity.NumberOfNodes;
SelSims.Settings.Seed=SimulationOptions.seed;

%% Save Simulation
% run DataExport.m
end
%run ShowMe.m