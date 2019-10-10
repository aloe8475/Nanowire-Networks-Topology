%Topology Analysis Simulation Script Start:


simType=lower(input('Choose Simulation: C - Continuous DC, P - 10 Pulse DC, T - Time Delay Analysis \n','s'));

numSims=100;
numTimes=100;
numNanowires=100;
inputVoltage=2;

switch simType
    case 'c'
        %% Continuous DC
        savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Single DC Pulse\';
        randseed=i;
        % contactn = randi(100,[1 2]);
        % timeDelay = i*0.05;
        biasType='DC';
        SelSims=runSimulation(contactn, [],randseed,biasType,numNanowires,inputVoltage);
    case 'p'
        %% 10 DC pulses: HAVE TO RUN THIS BEFORE TIME-DELAY (NEED TO USE SAME CONTACTN
        savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\';
        for i = 1:numSims
            randseed=i;
            contactn = randi(100,[1 2]);
            % timeDelay = i*0.05;
            biasType='DCandWait';
            SelSims{i}=runSimulation(contactn, [],randseed,biasType,numNanowires,inputVoltage);
        end
        save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
        
    case 't'
        %% 4 pulses with different Time-Delay bw pulse 3 and 4
        loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\10 Square Pulses\'; %load contact nodes from previous simulations

        savepath= 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
        load('Zdenka_100nw_DCandWait_100SimsOnly_2_Sec_2Electrodes_Vmax_2_09-Oct-2019.mat') %load contact nodes from previous simulations
        OldSims=SelSims; %save simulations from previous.
        clear SelSims
        for i = 1:numSims
            randseed=i;
            biasType='TimeDelay';
            contactn = [OldSims{i}.Electrodes.PosIndex]; %load contact nodes from previous simulation (case P)
            fprintf([num2str(i) '\n']);
            for j = 1:numTimes
                timeDelay = j*0.05;
                SelSims{i,j}=runSimulation(contactn, timeDelay,randseed,biasType,numNanowires,inputVoltage);
            end
        end
                save([savepath SelSims{1}.Settings.Model '_' num2str(SelSims{1}.NumberOfNodes) 'nw_' SelSims{1}.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims{1}.Settings.Time) '_Sec_' num2str(length(SelSims{1}.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims{1}.Settings.Vmax) '_' date '.mat'],'SelSims','-v7.3');
end
