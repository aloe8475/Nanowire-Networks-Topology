% Make sure you have added the folder 'asn', including its subdiretories
%close all;
close all;

numSims=100;
% pool=parpool(6);

for i = 1:numSims
%% Set the seed for PRNGs for reproducibility:
rng(42)
s = rng;
fprintf([num2str(i) '\n']);
%% Plot and analysis Output flags:
SimulationOptions{i}.takingsnapshots = true; % true \ false
SimulationOptions{i}.compilingMovie  = false; % true \ false 
SimulationOptions{i}.onlyGraphics    = false; % true \ false (no analysis is done and shown, only graphics (snapshots, movie) are generated).

%% Simulation general options:
SimulationOptions{i}.seed = s;    % save
SimulationOptions{i}.dt = 0.01;%1e-3;   % (sec)
SimulationOptions{i}.T  = 20;    % (sec) duration of simulation
SimulationOptions{i}.TimeVector = (SimulationOptions{i}.dt:SimulationOptions{i}.dt:SimulationOptions{i}.T)';
SimulationOptions{i}.NumberOfIterations = length(SimulationOptions{i}.TimeVector);  

%% Simulation recording options:
SimulationOptions{i}.ContactMode  = 'farthest';    % 'farthest' \ 'specifiedDistance' \ 'random' (the only one relevant for 'randAdjMat' (no spatial meaning)) \ 'preSet'
% SimulationOptions{i}.ContactNodes = [73, 51]; % only really required for preSet, other modes will overwrite this
%SimulationOptions{i}.isSource     = []; %only required for contacts with index >= 3


%% Generate Connectivity:
Connectivity{i}.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat' \ 'lattice'
switch Connectivity{i}.WhichMatrix
    case 'nanoWires'
        Connectivity{i}.filename = 'AdriantoZdenka100nw.mat';%'2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat'; %100nw
        Connectivity{i}.DataType = 'Adrian';
        %Connectivity{i}.filename  = '2016-09-08-153543_asn_nw_02048_nj_11469_seed_042_avl_28.00_disp_10.00.mat';
        %Connectivity{i}.filename = '2016-09-08-155044_asn_nw_00700_nj_14533_seed_042_avl_100.00_disp_10.00.mat'; %700nw
%         Connectivity{i}.filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
%         Connectivity{i}.DataType = 'Zdenka';
        %Connectivity{i}.filename = '2019-07-26-143006_asn_nw_05098_nj_11631_seed_042_avl_100.00_disp_10.00_gns_10.00_cdisp_1500.00.mat';
        %Connectivity{i}.filename = '2019-07-07-134601_asn_nw_00099_nj_00582_seed_501_avl_100.00_disp_10.00_gns_03.00_cdisp_100.00.mat';
    case 'randAdjMat'
        Connectivity{i}.NumberOfNodes = 30;
        Connectivity{i}.AverageDegree = 10;
    case 'Lattice'
        Connectivity{i}.sizex = 10;
        Connectivity{i}.sizey = 10;
end
Connectivity{i} = getConnectivity(Connectivity{i});

%% Choose  contacts:
if strcmp(SimulationOptions{i}.ContactMode, 'specifiedDistance')
    SimulationOptions{i}.BiProbeDistance = 500; % (um)1.5
end
SimulationOptions{i} = selectContacts(Connectivity{i}, SimulationOptions{i});


%% Initialize dynamic Components:
Components{i}.ComponentType       = 'atomicSwitch'; % 'atomicSwitch' \ 'memristor' \ 'resistor', \'tunnelSwitch' ,  \'quantCSwitch'
Components{i} = initializeComponents(Connectivity{i}.NumberOfEdges,Components{i});

%% Initialize Stimulus:
Stimulus{i}.BiasType              = 'AlonPulse';           % 'DC' \ 'AC' \ 'DCandWait' \ 'Ramp' \ 'ACsaw'
switch Stimulus{i}.BiasType
    case 'DC'
        Stimulus{i}.Amplitude = 0.2;%0.0448973; %0.045;  % (Volt)
    case {'AC', 'ACsaw'}
        Stimulus{i}.Frequency       = 1/5.0; % (Hz)
        Stimulus{i}.Amplitude       = 1.0;   % (Volt)
    case 'DCandWait'
        Stimulus{i}.OffTime      = 10.0; % SimulationOptions{i}.T/3; % (sec)
        Stimulus{i}.AmplitudeOn  = 0.2;                   % (Volt)
        Stimulus{i}.AmplitudeOff = 1e-3;                 % (Volt)
    case 'Ramp'
        Stimulus{i}.AmplitudeMin = 0;    % (Volt)
        Stimulus{i}.AmplitudeMax = 5;    % (Volt)
        
    case 'Triangle'
        Stimulus{i}.AmplitudeMin = 0;
        Stimulus{i}.AmplitudeMax = 3;    % (Volt)   
        
    case 'AlonPulse'
        Stimulus{i}.AmplitudeOn  = 1.5;
        Stimulus{i}.AmplitudeOff = 1e-3;
        Stimulus{i}.Period       = 3; %period of the short pulses
        Stimulus{i}.LongWait     = i*0.05; %Waiting time between the first set and second set of pulses
        Stimulus{i}.NumPulse1    = 3; %number of pulses before the long wait
        Stimulus{i}.NumPulse2    = 1; %number of pulses after the long wait
        
end
Stimulus{i} = getStimulus(Stimulus{i}, SimulationOptions{i});


%% Get Equations:
Equations{i} = getEquations(Connectivity{i},SimulationOptions{i}.ContactNodes,false);

%% Initialize snapshot time stamps:
if SimulationOptions{i}.takingsnapshots
    snapshotPeriod{i}   = 10*SimulationOptions{i}.dt; % (sec) make it a multiple integer of dt
    snapshotstep{i}     = ceil(snapshotPeriod{i} / SimulationOptions{i}.dt);
    snapshotsIdx{i}     = 1:snapshotstep{i}:SimulationOptions{i}.NumberOfIterations;
end

%% Simulate:
if SimulationOptions{i}.takingsnapshots
    [Output{i}, SimulationOptions{i}, snapshots{i}, SelSims{i}] = simulateNetwork(Connectivity{i}, Components{i}, Stimulus{i}, SimulationOptions{i}, snapshotsIdx{i}); % (Ohm)
else % this discards the snaphots
    [Output{i}, SimulationOptions{i}, snapshots{i}, SelSims{i}] = simulateNetwork(Connectivity{i}, Components{i}, Stimulus{i}, SimulationOptions{i}); % (Ohm)
end

%Convert Zdenka's structure to Adrian's Structure:
SelSims{i}=Convert_Zdenka_to_Adrian(SelSims{i},snapshots{i},SimulationOptions{i},Connectivity{i},Components{i},Stimulus{i});

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
% if ~SimulationOptions{i}.onlyGraphics
%     plotResults(Output{i}.networkResistance,Output{i}.networkCurrent,Stimulus{i});
% end


fprintf(' \n');
end
% Joel Hochstetter 15:52
save('VariableTime.mat','SelSims','Stimulus','Output');%% Graphics:
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