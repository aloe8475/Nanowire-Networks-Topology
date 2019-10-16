%Topology Analysis Simulation Script Start:

function callRunSimulation_PARALLEL(WorkerID)
%WorkerID is the parallel worker number which is assigned to this script

%% CHANGE SIMULATION TYPE HERE
% --------------------------------------------------------------------------------------------------------------------------
simType='t';%lower(input('Choose Simulation: C - Continuous DC, P - 10 Pulse DC, T - Time Delay Analysis, L - LDA \n','s'));
% --------------------------------------------------------------------------------------------------------------------------

%Initialise Variables
numSims      = 100;
numTimes     = 100;
numNanowires = 100;
inputVoltage =   2;
computer     = getenv('computername');


switch simType
    case 'c'
        %% Continuous DC
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Single DC Pulse\';
            case '' %if on linux
                savepath='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/Single DC Pulse/';
                %     case 'LAPTOP-S1BV3HR7'
                %         currentPath='D:\alon_\Research\PhD\CODE\Analysis';
                %case '' %--- Add other computer paths (e.g. Mike)
        end
        SimSettings.numSources=1;
        SimSettings.SimulationDuration=2;
        randseed = WorkerID;
        % contactn = randi(100,[1 2]);
        % timeDelay = i*0.05;
        biasType = 'DC';
        SelSims  = runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
        save([savepath 'SelSims_DC_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
    case 'p'
        %% 10 DC pulses: HAVE TO RUN THIS BEFORE TIME-DELAY (NEED TO USE SAME CONTACTN)
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\';
            case '' %if on linux
                savepath ='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
        end
        SimSettings.numSources=1;
        SimSettings.SimulationDuration=2;
        
        %         for i = WorkerID
        randseed   = WorkerID;
        contactn   = randi(100,[1 2]);
        % timeDelay = i*0.05;
        biasType   = 'DCandWait';
        SelSims{WorkerID} = runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
        %         end
        %         save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
        save([savepath 'SelSims_DCandWait_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
        fprintf('Data Saved');
%         %% Network Explore Analysis:
%         exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis/';
%         cd(exploreAnalysisPath);
%         exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/DC Pulse/';
%         [Explore{WorkerID}, threshold{WorkerID}, network{WorkerID}, Sim{WorkerID},analysis_type{WorkerID}]=Network_Explore_MultiSim_PARALLEL(numNanowires,WorkerID,biasType);
%         save([exploreSavePath 'ExploreAnalysis_DCandWait_' num2str(WorkerID) '.mat'],'Explore','threshold','network','Sim','analysis_type','-v7.3');
    case 't'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        %Initialise Paths
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
                loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
            case '' %if on linux
                savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
                loadpath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis/';
                exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/Time Delay Analysis/';
                
        end
        cd(savepath)
        if ~exist(['SelSims_TimeDelay_' num2str(WorkerID) '.mat'],'file') %if we haven't previously saved the file
            load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_1.5_09-Oct-2019.mat']) %load contact nodes from previous simulations
            SimSettings.numSources         = 1;
            SimSettings.SimulationDuration = 9;
            
            %         for i = WorkerID
            ElecPos(WorkerID,:) = [SelSims{WorkerID}.Electrodes.PosIndex]; %save electrodes from previous simulation
            %         end
            clear SelSims
            SelSims = cell(numSims,numTimes);
            %         for i = WorkerID
            randseed = WorkerID;
            biasType = 'TimeDelay'; % biasType for getStimulus function
            contactn = ElecPos(WorkerID,:); %load contact nodes from previous simulation (case P)
            fprintf([num2str(WorkerID) '\n']);
            for j = 1:numTimes
                timeDelay    = j*0.05;
                SelSims{WorkerID,j} = runSimulation(SimSettings,contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage);
            end
            %         end
            %         save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
            save([savepath 'SelSims_TimeDelay_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
            fprintf('Data Saved \n');
        else
            biasType = 'TimeDelay'; % biasType for getStimulus function
            fprintf('Data already exists... Loading \n');
            load([savepath 'SelSims_TimeDelay_' num2str(WorkerID) '.mat']);
            fprintf('Data Loaded \n');
        end
%         fprintf('Starting Explore Analysis \n');
%         %% Network Explore Analysis:
%         cd(exploreAnalysisPath);
%         [Explore{WorkerID}, threshold{WorkerID}, network{WorkerID}, Sim{WorkerID},analysis_type{WorkerID}]=Network_Explore_MultiSim_PARALLEL(numNanowires,WorkerID,biasType);
%         save([exploreSavePath 'ExploreAnalysis_TimeDelay_' num2str(WorkerID) '.mat'],'Explore','threshold','network','Sim','analysis_type','-v7.3');
%         
    case 'l'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4 + LDA Analysis Setup
        %         loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
        
        savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
        exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis/';
        exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/LDA Time Delay/';
        
        
        %         load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_2_09-Oct-2019.mat']) %load contact nodes from previous simulations
        
        SimSettings.numSources=2;
        SimSettings.SimulationDuration=9;
        %% NEED TO FIGURE OUT CONTACT NODE POSITIONS FOR LDA TOPOLOGY EXPLORATION - ALON 15/10/19
        contactn = randi(100,[1 4]); %Placement of 4 electrodes: [Source1, Drain1, Source2, Drain2]
        % 16/10/19 - ASK RUOMIN/JOEL/MIKE - NEED TO FIGURE OUT HOW TO TEST THE CONTACT NODES AND SELECT THEM
        % BASED ON THE DIFFERENT PARAMETERS: LowCOMM/LowCOMM, Low/Med, Low/High, Med/Low, Med/Med, Med/High, High/Low, High/Med, High/High.
        % OR can use Early/Mid/Late times from Experiment 2.
        
        clear SelSims
        SelSims=cell(numSims,numTimes);
        %         for i = WorkerID
        randseed=WorkerID;
        biasType{1}='TimeDelay1'; %different biasTypes for getStimulus function
        biasType{2}='TimeDelay2';
        tempcontact = contactn(WorkerID,:); %load contact nodes
        fprintf([num2str(WorkerID) '\n']);
        for j = 1:numTimes
            timeDelay = j*0.05;
            SelSims{WorkerID,j}=runSimulation(SimSettings,tempcontact, timeDelay,randseed,biasType,numNanowires,inputVoltage);
        end
        clear tempcontact
        %         end
        %         save([savepath SelSims{1,1}.Settings.Model '_' num2str(SelSims{1,1}.NumberOfNodes) 'nw_' SelSims{1,1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1,1}.Settings.Time) '_Sec_' num2str(length(SelSims{1,1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1,1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
        save([savepath 'SelSims_TimeDelayLDA_' num2str(WorkerID) '.mat'],'SelSims','-v7.3');
        fprintf('Data Saved');
%         %% Network Explore Analysis:
%         cd(exploreAnalysisPath);
%         [Explore{WorkerID}, threshold{WorkerID}, network{WorkerID}, Sim{WorkerID},analysis_type{WorkerID}]=Network_Explore_MultiSim_PARALLEL(numNanowires,WorkerID,biasType);
%         save([exploreSavePath 'ExploreAnalysis_TimeDelayLDA_' num2str(WorkerID) '.mat'],'Explore','threshold','network','Sim','analysis_type','-v7.3');
end

end