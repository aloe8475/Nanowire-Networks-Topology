% Make sure you have added the folder 'asn', including its subdiretories
%close all;
close all;

numSims=100;
% pool=parpool(6);
function runSimulation(contactn, timeDelay)
%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;
fprintf([num2str '\n']);
%% Plot and analysis Output flags:
SimulationOptions.takingsnapshots = true; % true \ false
SimulationOptions.compilingMovie  = false; % true \ false 
SimulationOptions.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions.seed = s;    % save
SimulationOptions.dt = 0.01;%1e-3;   % (sec)
SimulationOptions.T  = 20;    % (sec) duration of simulation
SimulationOptions.TimeVector = (SimulationOptions.dt:SimulationOptions.dt:SimulationOptions.T)';
SimulationOptions.NumberOfIterations = length(SimulationOptions.TimeVector);  

%% Simulation recording options:
SimulationOptions.ContactMode  = 'preSet';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
SimulationOptions.ContactNodes = contactn; % only really required for preSet, other modes will overwrite this
%SimulationOptions.isSource     = []; %only required for contacts with index >= 3


%% Generate Connectivity:
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat' \ 'lattice'
switch Connectivity.WhichMatrix
    case 'nanoWires'
        Connectivity.filename = 'AdriantoZdenka100nw_simulation1.mat';%'2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat'; %100nw
        Connectivity.DataType = 'Adrian';
        %Connectivity.filename  = '2016-09-08-153543_asn_nw_02048_nj_11469_seed_042_avl_28.00_disp_10.00.mat';
        %Connectivity.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat'; %700nw
%         Connectivity.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
%         Connectivity.DataType = 'Zdenka';
        %Connectivity.filename = '2019-07-26-143006_asn_nw_05098_nj_11631_seed_042_avl_100.00_disp_10.00_gns_10.00_cdisp_1500.00.mat';
        %Connectivity.filename = '2019-07-07-134601_asn_nw_00099_nj_00582_seed_501_avl_100.00_disp_10.00_gns_03.00_cdisp_100.00.mat';
    case 'randAdjMat'
        Connectivity.NumberOfNodes = 30;
        Connectivity.AverageDegree = 10;
    case 'Lattice'
        Connectivity.sizex = 10;
        Connectivity.sizey = 10;
end
Connectivity = getConnectivity(Connectivity);

%% Choose  contacts:
if strcmp(SimulationOptions.ContactMode, 'specifiedDistance')
    SimulationOptions.BiProbeDistance = 500; % (um)1.5
end
SimulationOptions = selectContacts(Connectivity, SimulationOptions);


%% Initialize dynamic Components:
Components.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor', \'tunnelSwitch' ,  \'quantCSwitch'
Components = initializeComponents(Connectivity.NumberOfEdges,Components);

%% Initialize Stimulus:
Stimulus.BiasType              = 'AlonPulse';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'ACsaw'
switch Stimulus.BiasType
    case 'DC'
        Stimulus.Amplitude = 0.2;%0.0448973; %0.045;  % (Volt)
    case {'AC', 'ACsaw'}
        Stimulus.Frequency       = 1/5.0; % (Hz)
        Stimulus.Amplitude       = 1.0;   % (Volt)
    case 'DCandWait'
        Stimulus.OffTime      = 10.0; % SimulationOptions.T/3; % (sec)
        Stimulus.AmplitudeOn  = 0.2;                   % (Volt)
        Stimulus.AmplitudeOff = 1e-3;                 % (Volt)
    case 'Ramp'
        Stimulus.AmplitudeMin = 0;    % (Volt)
        Stimulus.AmplitudeMax = 5;    % (Volt)
        
    case 'Triangle'
        Stimulus.AmplitudeMin = 0;
        Stimulus.AmplitudeMax = 3;    % (Volt)   
        
    case 'AlonPulse'
        Stimulus.AmplitudeOn  = 3;
        Stimulus.AmplitudeOff = 1e-3;
        Stimulus.Period       = 3; %period of the short pulses
        Stimulus.LongWait     = timeDelay; %Waiting time between the first set and second set of pulses
        Stimulus.NumPulse1    = 3; %number of pulses before the long wait
        Stimulus.NumPulse2    = 1; %number of pulses after the long wait
        
end
Stimulus = getStimulus(Stimulus, SimulationOptions);


%% Get Equations:
Equations = getEquations(Connectivity,SimulationOptions.ContactNodes,false);

%% Initialize snapshot time stamps:
if SimulationOptions.takingsnapshots
    snapshotPeriod   = 10*SimulationOptions.dt; % (sec) make it a multiple integer of dt
    snapshotstep     = ceil(snapshotPeriod / SimulationOptions.dt);
    snapshotsIdx     = 1:snapshotstep:SimulationOptions.NumberOfIterations;
end

%% Simulate:
if SimulationOptions.takingsnapshots
    [Output, SimulationOptions, snapshots, SelSims] = simulateNetworkMulti(Connectivity, Components, Stimulus, SimulationOptions, snapshotsIdx); % (Ohm)
else % this discards the snaphots
    [Output, SimulationOptions, snapshots, SelSims] = simulateNetworkMulti(Connectivity, Components, Stimulus, SimulationOptions); % (Ohm)
end

%Convert Zdenka's structure to Adrian's Structure:
SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions,Connectivity,Components,Stimulus);

% Connectivity Edge Position for Network with Simulation
 %NEED TO CHANGE EDGE POISITION TO NODE POSITION 
 % (i.e. inverse this:)
 %
% xi=triu(LayoutSim.CX);
% xi=xi(adj_matrix~=0);
% xi=xi(xi~=0)';
% 
% yi=triu(LayoutSim.CY);
% yi=yi(adj_matrix~=0);
% yi=yi(yi~=0)';

% LayoutSim.CX=xi;
% LayoutSim.CY=yi; 

%% Analysis and plot results:
% if ~SimulationOptions.onlyGraphics
%     plotResults(Output.networkResistance,Output.networkCurrent,Stimulus);
% end


fprintf(' \n');

% Joel Hochstetter 15:52
save('VariableTime.mat','SelSims','Stimulus','-v7.3');%% Graphics:
end 
% if SimulationOptions.takingsnapshots
%     % What to plot:
%     whatToPlot = struct(...
%                         'Nanowires',    true, ...
%                         'Contacts',     true, ...
%                         'OnOrOff',      true, ...
%                         'Dissipation',  false, ...
%                         'Lambda',       true, ... #can either plot lambda or dissipation or Vdrop
%                         'Currents',     true, ...
%                         'Voltages',     true, ...
%                         'Labels',       true, ...
%                         'VDrop',        true, ... 
%                         'GraphRep',     true ...
%                         );
%     
%                     
%     % Uniform scales:
%     axesLimits.DissipationCbar = [0,5]; % (1pW*10^0 : 1pW*10^5)
%     axesLimits.CurrentArrowScaling = 1.1;
%     axesLimits.LambdaCbar = [0,0.15];
%     switch Stimulus.BiasType
%         case 'DC'
%             if Stimulus.Signal(1) > 0
%                 axesLimits.VoltageCbar = [0,Stimulus.Signal(1)]; % (V)
%             else
%                 axesLimits.VoltageCbar = [Stimulus.Signal(1),0]; % (V)
%             end
%         case {'AC' , 'DCandWait', 'Ramp', 'Triangle', 'ACsaw'}
%             axesLimits.VoltageCbar = [min(Stimulus.Signal),max(Stimulus.Signal)]; % (V)
%     end
% 	
%     
%     % Just for fun - extract a specific frame (from the middle):
%     %snapshotToFigure(snapshots{floor(length(snapshots)/2)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%     snapshotToFigure(snapshots{end},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%     set(gcf, 'visible','on')
% 
%     % Compile whole movie:
%     if SimulationOptions.compilingMovie 
%         fprintf('\nCompiling movie...\n');
%         
%         % Only Windows and MacOS >  10.7
%         %v = VideoWriter('networkMovie.mp4','MPEG-4');
%         % All platformse.g. for DC this is just the minimum |lambda-lambda_c|/|V-V_set| for all the switches)
%         v = VideoWriter('networkMovie','Motion JPEG AVI');
% 
%         v.FrameRate = floor(1/snapshotPeriod/10);
%         
%         v.Quality = 100;
%         open(v);
%         for i = 1 : length(snapshots)
%             progressBar(i,length(snapshots));
%             frameFig = snapshotToFigure(snapshots,SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%             writeVideo(v,getframe(frameFig));
%             close(frameFig);
%         end
%         close(v);
%         fprintf('\nDone.\n');
%     end
% end
% 
% fprintf('\n');