% 20/04/19
% Loading Function for Adrian's Code v0.1

function [network, sim_loaded, explore_network, numNetworks] = Load_Adrian_Code()
%% Load Network Data:
%choose network to load
computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        dataPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Adrian Networks';
    case ''
        dataPath='/suphys/aloe8475/Documents/CODE/Adrian Code/NETWORK_sims_2/Saved Networks';
    case 'LAPTOP-S1BV3HR7'
        dataPath='D:\alon_\Research\PhD\CODE\Adrian Code\NETWORK_sims_2\Saved Networks';
end
cd(dataPath)
waitfor(msgbox('Select the Network saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
count=1;
f{count}=fullfile(PathName,FileName);
load(f{count});

% Define network variable
network(count)=SelNet;
clear SelNet
% Delete default simulation data:

numNetworks=1;
i=1;
more_networks=input('Do you want to load another Network? \n','s');
if more_networks=='y'
    while i==1
        %Update network count
        count=count+1;
        %prompt for the network .mat file with button
        promptMessage = sprintf('Select the Network saved data');
        button = questdlg(promptMessage, 'Load Network Data', 'OK','Cancel','OK');
        if strcmp(button, 'Cancel')
            numNetworks=count-1; %if they change their mind and cancel, skip to simulations
            break;
        elseif strcmp(button,'OK')
            % if they want a second network, prompt for second network data
            [FileName2,PathName2] = uigetfile('*.mat','Select the Network saved data');
            f{count}=fullfile(PathName2,FileName2);
            if count ~= 1 & sum(~cellfun('isempty',strfind(f,f{count})))>1 %if the current network string occurs more than once
                %ask them to select a different network
                waitfor(msgbox('Error: Cannot load same network twice. Please select a different Network'));
                count=count-1;%reduce our counter because we go back to the start of the while loop. 
                error_state=1;
            else
                load(f{count});
                % Define network variable
                network(count)=SelNet;
                % save a var with number of networks
                numNetworks=count;
                i=i+1; %can also use break
                error_state=0;
            end
        end
        if error_state~=1 %we don't want to ask to load again if we got an error
        more_networks=input('Do you want to load another Network? \n','s');
        else
            more_networks='y';
        end 
        if more_networks=='y'
            i=1;
        end 
    end
end

%% Load Simulation Data
sims_load=input('Load Simulation Data too? y or n? \n','s');
if sims_load=='y'
    % Delete default simulation data:
    network(1).Simulations = [];
    if length(network)>1
             % Delete default simulation data as well:
                network(2).Simulations = [];
    end 
    explore_network=lower(input('Do you want Training + Testing Data, or Just Explore Data? - T or E \n','s'));
    %Select which Simulation to Load
    sim_loaded=1;
    if explore_network=='t'
        cd('Simulations Only\Training');
        waitfor(msgbox('Select the Training Simulation saved data'))
        [FileName,PathName] = uigetfile('*.mat','Select the Training saved data');
        f=fullfile(PathName,FileName);
        load(f);
        
        %Add simulation data to network struct:
        if numNetworks==2
            training_network=input('Which Network was the Training simulation from? 1 or 2 \n');
        end
        temp1=SelSims;
        
        clear SelSims
        cd('..\Testing');
        waitfor(msgbox('Select the Testing Simulation saved data'))
        [FileNameTest,PathNameTest] = uigetfile('*.mat','Select the Testing saved data');
        f_test=fullfile(PathNameTest,FileNameTest);
        load(f_test);
        
        %Add simulation data to network struct:
        if numNetworks==2
            testing_network=input('Which Network was the Testing simulation from? 1 or 2 \n');
        end
        if exist('SelSims','var') == 1
        temp2 = SelSims;%save simulations from test network if it hasn't been changed
        else 
        temp2 = network.Simulations;%save simulations from test network if it has been changed
        end 
        
        if numNetworks==1 %if both training and testing are from same networks, just combine the simulations
            network.Simulations=[temp1 temp2];
            for j = 1:length(temp1)
            network.Simulations{j}.Type='Training Simulation'; %label the training sim
            network.numTrainingSims=length(temp1);
            end 
            for i=(length(temp1)+1):length(network.Simulations)
            network.Simulations{i}.Type='Testing Simulation';
            network.numTestingSims=length(network.Simulations)-length(temp1);
            end 
        elseif numNetworks==2 %if the two networks are different, save the simulations in the different networks
            network(training_network).Simulations=temp1;
            network(testing_network).Simulations=temp2;
            network(training_network).Simulations{1}.Type='Training Simulation'; %label which is the training sim
            network(training_network).Type='Training Network';
            clear SelSims
        end
    else
        cd('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Adrian Only\Explore'); %change for all computers
        waitfor(msgbox('Select the Explore Simulation saved data'))
        [FileName,PathName] = uigetfile('*.mat','Select the Explore saved data');
        f_explore=fullfile(PathName,FileName);
        load(f_explore);
        if exist('SelSims','var') == 1
        temp= SelSims;%save simulations from test network if it hasn't been changed
        else 
        temp = network.Simulations;%save simulations from test network if it has been changed
        end 
        network.Simulations=temp;
        for i=1:length(network.Simulations)
        network.Simulations{i}.Type='Explore Simulation';
        end 
        %Need to change to only allow explore of 1 network
        clear SelSims;
    end
else
    sim_loaded=0;
    explore_network='e';
end
fprintf('Network and Simulations Loaded \n');
fprintf('\n--------------------------------- \n\n');
end

