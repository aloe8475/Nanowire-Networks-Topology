% 20/04/19
% Loading Function for Adrian's Code v0.1

function [network, sim_loaded] = Load_Adrian_Code()
%% Load Network Data:
%choose network to load
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks')
waitfor(msgbox('Select the Network saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
f=fullfile(PathName,FileName);
load(f);

% Define network variable
network=SelNet;
clear SelNet
% Delete default simulation data:
network.Simulations = [];

%% Load Simulation Data
sims_load=input('Load Simulation Data too? y or n? \n','s');
if sims_load=='y'
    %Select which Simulation to Load
    sim_loaded=1;
    cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only');
    waitfor(msgbox('Select the Simulation saved data'))
    [FileName,PathName] = uigetfile('*.mat','Select the  saved data');
    f=fullfile(PathName,FileName);
    load(f);
    
    %Add simulation data to network struct:
    network.Simulations = SelSims;
    clear SelSims
else
        sim_loaded=0;
end 
end

