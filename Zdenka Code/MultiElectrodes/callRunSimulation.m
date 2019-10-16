%Topology Analysis Simulation Script Start:


simType=lower(input('Choose Simulation: C - Continuous DC, P - 10 Pulse DC, T - Time Delay Analysis, L - LDA \n','s'));

numSims=1;
numTimes=100;
numNanowires=100;
inputVoltage=10;
computer=getenv('computername');

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
        randseed=i;
        % contactn = randi(100,[1 2]);
        % timeDelay = i*0.05;
        biasType='DC';
        SelSims=runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
    case 'p'
        %% 10 DC pulses: HAVE TO RUN THIS BEFORE TIME-DELAY (NEED TO USE SAME CONTACTN
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\';
            case '' %if on linux
                savepath='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
        end
        SimSettings.numSources=1;
        SimSettings.SimulationDuration=2;
        
        for i = 1:numSims
            randseed=i;
            contactn = randi(100,[1 2]);
            % timeDelay = i*0.05;
            biasType='DCandWait';
            SelSims{i}=runSimulation(SimSettings,contactn, [],randseed,biasType,numNanowires,inputVoltage);
        end
        save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
                fprintf('Data Saved');

    case 't'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
                loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
            case '' %if on linux
                savepath='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
                loadpath='/suphys/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/10 Square Pulses/';
                
        end
        load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_2_09-Oct-2019.mat']) %load contact nodes from previous simulations
        SimSettings.numSources=1;
        SimSettings.SimulationDuration=9;
        
        for i = 1:numSims
            ElecPos(i,:)=[SelSims{i}.Electrodes.PosIndex]; %save simulations from previous.
        end
        clear SelSims
        SelSims=cell(numSims,numTimes);
        for i = 1:numSims
            randseed=i;
            biasType='TimeDelay';
            contactn = ElecPos(i,:); %load contact nodes from previous simulation (case P)
            fprintf([num2str(i) '\n']);
            for j = 1:numTimes
                timeDelay = j*0.05;
                SelSims{i,j}=runSimulation(SimSettings,contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage);
            end
        end
        save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
        fprintf('Data Saved');
    case 'l'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        %         loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations
        
        savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
        %         load([loadpath 'Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_2_09-Oct-2019.mat']) %load contact nodes from previous simulations
        
        SimSettings.numSources=2;
        SimSettings.SimulationDuration=9;
        contactn = randi(100,[1 4]); %Placement of 4 electrodes: [Source1, Drain1, Source2, Drain2]
        clear SelSims
        SelSims=cell(numSims,numTimes);
        for i = 1:numSims
            randseed=i;
            biasType{1}='TimeDelay1';
            biasType{2}='TimeDelay2';
            tempcontact = contactn(i,:); %load contact nodes from previous simulation (case P)
            fprintf([num2str(i) '\n']);
            for j = 1:numTimes
                timeDelay = j*0.05;
                SelSims{i,j}=runSimulation(SimSettings,tempcontact, timeDelay,randseed,biasType,numNanowires,inputVoltage);
            end
            clear tempcontact
        end
        save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
                fprintf('Data Saved');

        
end
