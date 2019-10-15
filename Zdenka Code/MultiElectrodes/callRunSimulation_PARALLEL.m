%Topology Analysis Simulation Script Start:

function callRunSimulation_PARALLEL(WorkerID)
simType='t';%lower(input('Choose Simulation: C - Continuous DC, P - 10 Pulse DC, T - Time Delay Analysis, L - LDA \n','s'));

numSims      = 100;
numTimes     = 100;
numNanowires = 100;
inputVoltage =   2;
computer     = getenv('computername');

switch simType
    case 'c'
        %% Continuous DC
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
    case 'p'
        %% 10 DC pulses: HAVE TO RUN THIS BEFORE TIME-DELAY (NEED TO USE SAME CONTACTN
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
        %% Network Explore Analysis:
        exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis';
        cd(exploreAnalysisPath);
        exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/DC Pulse/';
        Explore{WorkerID}=Network_Explore_MultiSim_PARALLEL(numNanowires,SelSims{WorkerID}); 
        save([exploreSavePath 'ExploreAnalysis_TimeDelay_' num2str(WorkerID) '.mat'],'Explore','-v7.3');
    case 't'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
                loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
            case '' %if on linux
                savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
                loadpath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                
        end
        load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_1.5_09-Oct-2019.mat']) %load contact nodes from previous simulations
        SimSettings.numSources         = 1;
        SimSettings.SimulationDuration = 9;
        
        %         for i = WorkerID
        ElecPos(WorkerID,:) = [SelSims{WorkerID}.Electrodes.PosIndex]; %save simulations from previous.
        %         end
        clear SelSims
        SelSims = cell(numSims,numTimes);
        %         for i = WorkerID
        randseed = WorkerID;
        biasType = 'TimeDelay';
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
        fprintf('Starting Explore Analysis');
        %% Network Explore Analysis:
        exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis';
        cd(exploreAnalysisPath);
        exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/Time Delay Analysis/';
        Explore{WorkerID}=Network_Explore_MultiSim_PARALLEL(numNanowires,SelSims{WorkerID,:}); 
        save([exploreSavePath 'ExploreAnalysis_TimeDelay_' num2str(WorkerID) '.mat'],'Explore','-v7.3');
        
    case 'l'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        %         loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
        
        savepath ='/headnode2/aloe8475/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
        %         load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_2_09-Oct-2019.mat']) %load contact nodes from previous simulations
        
        SimSettings.numSources=2;
        SimSettings.SimulationDuration=9;
        contactn = randi(100,[1 4]); %Placement of 4 electrodes: [Source1, Drain1, Source2, Drain2]
        clear SelSims
        SelSims=cell(numSims,numTimes);
        %         for i = WorkerID
        randseed=WorkerID;
        biasType{1}='TimeDelay1';
        biasType{2}='TimeDelay2';
        tempcontact = contactn(WorkerID,:); %load contact nodes from previous simulation (case P)
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
         %% Network Explore Analysis:
        exploreAnalysisPath='/headnode2/aloe8475/CODE/Analysis';
        cd(exploreAnalysisPath);
        exploreSavePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/LDA Time Delay/';
        Explore{WorkerID}=Network_Explore_MultiSim_PARALLEL(numNanowires,SelSims{WorkerID,:}); 
        save([exploreSavePath 'ExploreAnalysis_TimeDelay_' num2str(WorkerID) '.mat'],'Explore','-v7.3');
end

end
