% Make sure you have added the folder 'asn', including its subdiretories
close all;
clearvars -except Nodes testcurrent testconduct contactn recurrent

tic
%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;

numSims=1;
for simNum=1:numSims

    
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
SimulationOptions.ContactMode     = 'farthest';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
% FOR MULTIPLE ELECTRODES NEED TO DEFINE BELOW AND UNCOMMENT BOTH LINES
% SimulationOptions.electrodes      = [73,30];
if strcmp(SimulationOptions.ContactMode,'preSet')
SimulationOptions.numOfElectrodes = length(SimulationOptions.electrodes);
else
SimulationOptions.numOfElectrodes = 2;
end 


%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
switch Connectivity.WhichMatrix
    case 'nanoWires'
          numNanowires=100c; %Change size to load different networks 
          savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\';
          loadpath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Zdenka Networks\';
          Connectivity.filename = [loadpath 'AdriantoZdenka' num2str(numNanowires) 'nw_simulation1.mat'];
%         Connectivity.filename = 'NetworkData/2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
        %Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
end
Connectivity = getConnectivityMulti(Connectivity);

%% Choose  contacts:


if ~strcmp(SimulationOptions.ContactMode,'preSet')
    if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 1500; % (um)
    end
    SimulationOptions = selectContacts(Connectivity, SimulationOptions);
end 

%% Initialize dynamic components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor'
Components = initializeComponentsMulti(Connectivity.NumberOfEdges,Components);

%% Initialize stimulus:
Signals = cell(SimulationOptions.numOfElectrodes,1);

Stimulus1.BiasType       = 'DC';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'AlonPulse'
% Stimulus1.OnTime         = 0.0;
% Stimulus1.OffTime        = SimulationOptions.T/10; %ten pulses
Stimulus1.AmplitudeOn    = 1;
Stimulus1.AmplitudeOff   = 0.005;
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
if SimulationOptions.takingSnapshots
[Output, SimulationOptions, snapshots, SelSims{simNum},testcurrent,testconduct] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions, snapshotsIdx); % (Ohm)
else % this discards the snaphots
[Output, SimulationOptions, snapshots, SelSims{simNum},testcurrent,testconduct] = simulateNetworkMulti(Connectivity, Components, Signals, SimulationOptions); % (Ohm)
end


%Convert Zdenka's structure to Adrian's Structure:
SelSims{simNum}=Convert_Zdenka_to_Adrian(SelSims{simNum},snapshots,SimulationOptions,Connectivity,Components,Stimulus);
SelSims{simNum}.Settings.SigType = Stimulus{1}.BiasType;
fprintf('\n')
end 
%% Save Simulation
save([savepath SelSims{simNum}.Settings.Model '_' num2str(Connectivity.NumberOfNodes) 'nw_' SelSims{simNum}.Settings.SigType '_' num2str(length(SelSims{simNum})) 'SimsOnly_' num2str(SelSims{simNum}.Settings.Time) '_Sec_' num2str(length(SelSims{simNum}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{simNum}.Settings.Vmax) '_' date],'SelSims','-v7.3');
% run DataExport.m
toc
%run ShowMe.m 