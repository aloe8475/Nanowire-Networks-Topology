%Topology Analysis Simulation Script Start:

function callRunSimulation_PARALLEL(WorkerID)
%WorkerID is the parallel worker number which is assigned to this script

%% CHANGE SIMULATION TYPE HERE
% -------------------------------------------------------------------------------------------------------------------------
simType='c';%lower(input('Choose Simulation: C - Continuous DC, P - 10 Pulse DC, T - Time Delay Analysis, L - LDA \n','s'));
pathLengths='d'; %s = same, d = different;
% --------------------------------------------------------------------------------------------------------------------------

%Initialise Variables
numSims      = 100;
numTimes     = 200;
numNanowires = 100;
inputVoltage =   1;
computer     = getenv('computername');
switch computer
    case 'W4PT80T2'
        loadpairing='C:\Users\aloe8475\Dropbox (Sydney Uni)\Data\ASN_simulation\Python\ASN\connectivity_data\';
        
    case '' %if on linux
        loadpairing='/headnode2/aloe8475/CODE/Analysis/data/';
end

if pathLengths=='s'
    load([loadpairing 'ElecPosPathLength.mat']); %load electrode pairing positions.
    elecPos=ElecPos+1;
    clear ElecPos
else
    load([loadpairing 'ElecPos.mat']); %load electrode pairing positions.
    % load([loadpairing 'ElecPosPathLength.mat']);%- IF WE WANT ALL OF SAME PATHLENGTH
    %Modify some electrode positions that were the same as others:
    elecPos(1,:)=[5,18];
    elecPos(81,:) =[23, 81];
    elecPos(82,:) = [100, 82];
end
switch simType
    case 'c'
        %% Continuous DC
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Single DC Pulse\';
            case '' %if on linux
                %                                 savepath='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/Single DC Pulse/';
                savepath='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Single DC Pulse/';
                %     case 'LAPTOP-S1BV3HR7'
                %         currentPath='D:\alon_\Research\PhD\CODE\Analysis';
                %case '' %--- Add other computer paths (e.g. Mike)
        end
       
        SimSettings.numSources=1;
        SimSettings.SimulationDuration=2;
        randseed = WorkerID*2;
        contactn = elecPos(WorkerID,:);
        % timeDelay = i*0.05;
        biasType{1} = 'DC';
        SelSims  = runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
        save([savepath 'SelSims_DC_' pathLengths '_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
    case 'p'
        %% 10 DC pulses: HAVE TO RUN THIS BEFORE TIME-DELAY (NEED TO USE SAME CONTACTN)
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\';
                exploreAnalysisPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\';
                exploreSavePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\DC Pulse\';
                
            case '' %if on linux
                %                 savepath ='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis/';
                exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/DC Pulse/';
                
                
        end
        cd(savepath)
        if ~exist(['SelSims_DCandWait_' pathLengths '_' num2str(WorkerID) '.mat'],'file') %if we haven't previously saved the file
            
            SimSettings.numSources=1;
            SimSettings.SimulationDuration=2;
            
            %         for i = WorkerID
            randseed   = WorkerID;
            contactn   = elecPos(WorkerID,:);
            % timeDelay = i*0.05;
            biasType{1}   = 'DCandWait';
            SelSims{WorkerID} = runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
            %         end
            fprintf('Starting Data Analysis \n');
            %         save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
            save([savepath 'SelSims_DCandWait_' pathLengths '_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
            fprintf('Data Saved \n');
        else
            biasType{1} = 'DCandWait'; % biasType for getStimulus function
            fprintf('Data already exists... Loading \n');
            load([savepath 'SelSims_DCandWait_' pathLengths '_' num2str(WorkerID) '.mat']);
            fprintf('Data Loaded \n');
        end
        fprintf('Starting DC Pulse Explore Analysis \n');
        
        %% Network Explore Analysis:
        cd(exploreAnalysisPath);
        [Explore{WorkerID}, threshold{WorkerID}, network{WorkerID}, Sim{WorkerID},analysis_type{WorkerID}]=Network_Explore_MultiSim_PARALLEL(numNanowires,WorkerID,biasType,SelSims{WorkerID});
        save([exploreSavePath 'ExploreAnalysis_DCandWait_' pathLengths '_' num2str(WorkerID) '.mat'],'Explore','threshold','network','Sim','analysis_type','-v7.3');
        fprintf('Explore Saved');
        
    case 't'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
                loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
                exploreAnalysisPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\';
                exploreSavePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Time Delay Analysis\';
                
            case '' %if on linux
                linux=1; %% 0 = LVM, 1 = Cluster
                if linux==1
                    savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
                    loadpath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                    exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis/';
                    exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/Time Delay Analysis/';
                else
                    savepath ='/import/silo2/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
                    loadpath ='/import/silo2/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                    exploreAnalysisPath='/import/silo2/aloe8475/Documents/CODE/Analysis/';
                    exploreSavePath='/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/Time Delay Analysis/';
                end
        end
        cd(savepath)
        if ~exist(['SelSims_TimeDelay_' pathLengths '_' num2str(WorkerID) '.mat'],'file') %if we haven't previously saved the file
            SimSettings.numSources         = 1;
            SimSettings.SimulationDuration = (numTimes/10)+2;
            SelSims = cell(numSims,numTimes);
            %         for i = WorkerID
            randseed = WorkerID;
            biasType{1} = 'TimeDelay'; % biasType for getStimulus function
            contactn = elecPos(WorkerID,:); %ElecPos(WorkerID,:); %load contact nodes from previous simulation (case P)
            fprintf([num2str(WorkerID) '\n']);
            for j = 1:numTimes
                timeDelay    = j*0.05;
                fprintf(['Running TimeDelay ' num2str(timeDelay) '...'])
                SelSims{WorkerID,j} = runSimulation(SimSettings,contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage);
                fprintf(['\n'])
            end
            %         end
            %         save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
            save([savepath 'SelSims_TimeDelay_' pathLengths '_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
            % %             fprintf('Data Saved \n');
        else
            biasType{1} = 'TimeDelay'; % biasType for getStimulus function
            fprintf('Data already exists... Loading \n');
            load([savepath 'SelSims_TimeDelay_' pathLengths '_' num2str(WorkerID) '.mat']);
            fprintf('Data Loaded \n');
        end
        fprintf('Starting Time Delay Explore Analysis \n');
        %% Network Explore Analysis:
        cd(exploreAnalysisPath);
        %         if ~exist(['ExploreAnalysis_TimeDelay_' num2str(WorkerID) '.mat'],'file') %if we haven't previously saved the file
        allSims={SelSims{WorkerID,:}};
        [ExploreTemp{WorkerID}, threshold{WorkerID}, network{WorkerID}, Sim{WorkerID},analysis_type{WorkerID}]=Network_Explore_MultiSim_PARALLEL(numNanowires,WorkerID,biasType,allSims);
        for i = 1:length(ExploreTemp{WorkerID})
            Explore{WorkerID,i}=ExploreTemp{WorkerID}{i}; %save each time point & each workerID together
            %             a = Sim{WorkerID}{i}.Data;
            %             Sim{WorkerID}{i}.Data=rmfield(a,{'Rmat','IDrain1','VDrain1','WireCurrents','WireVoltages','ElectrodeCurrents','JunctionVoltages','JunctionResistance','JunctionCurrents','JunctionRmat','JunctionVoltage'});
            %             Sim{WorkerID}{i}=rmfield(Sim{WorkerID}{i},{'Gmat','Time','SelLayout'});
        end
        save([exploreSavePath 'ExploreAnalysis_TimeDelay_' pathLengths '_' num2str(WorkerID) '.mat'],'Explore','threshold','analysis_type','-v7.3');
        fprintf('Explore Saved');
        %         else
        %         fprintf('Explore Data already exists \n');
        %         end
        
end

end
